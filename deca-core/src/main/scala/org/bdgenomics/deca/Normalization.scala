/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.bdgenomics.deca

import breeze.linalg.sum
import breeze.stats.mean
import breeze.numerics.sqrt
import breeze.stats.meanAndVariance
import org.apache.spark.SparkContext
import org.apache.spark.mllib.linalg.SingularValueDecomposition
import org.apache.spark.mllib.linalg.distributed.{ IndexedRow, IndexedRowMatrix }
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.deca.Timers._
import org.bdgenomics.deca.coverage.ReadDepthMatrix
import org.bdgenomics.deca.util.MLibUtils
import org.bdgenomics.utils.misc.Logging

import scala.annotation.tailrec
import scala.util.control.Breaks._

object Normalization extends Serializable with Logging {

  def filterColumns(matrix: IndexedRowMatrix, targets: Array[ReferenceRegion],
                    minMean: Double = Double.MinValue,
                    maxMean: Double = Double.MaxValue,
                    minSD: Double = 0,
                    maxSD: Double = Double.MaxValue): (IndexedRowMatrix, Array[ReferenceRegion]) = {
    matrix.rows.cache()
    val colStats = matrix.toRowMatrix.computeColumnSummaryStatistics
    val colMeans = MLibUtils.mllibVectorToDenseBreeze(colStats.mean)
    val colSDs = sqrt(MLibUtils.mllibVectorToDenseBreeze(colStats.variance))

    val toKeep = (colMeans :>= minMean) :& (colMeans :<= maxMean) :& (colSDs :>= minSD) :& (colSDs :<= maxSD)
    val newTargets = targets.zipWithIndex.collect { case (target, index) if toKeep(index) => target }

    val toKeepBroadcast = SparkContext.getOrCreate().broadcast(toKeep)
    val filteredMatrix = new IndexedRowMatrix(matrix.rows.map(row => {
      val vector = MLibUtils.mllibVectorToDenseBreeze(row.vector)
      IndexedRow(row.index, MLibUtils.breezeVectorToMLlib(vector(toKeepBroadcast.value)))
    }))

    (filteredMatrix, newTargets)
  }

  def meanCenterColumns(matrix: IndexedRowMatrix): IndexedRowMatrix = {
    matrix.rows.cache()
    val colStats = matrix.toRowMatrix.computeColumnSummaryStatistics
    new IndexedRowMatrix(matrix.rows.map(row => {
      IndexedRow(
        row.index,
        MLibUtils.breezeVectorToMLlib(MLibUtils.mllibVectorToDenseBreeze(row.vector) - MLibUtils.mllibVectorToDenseBreeze(colStats.mean)))
    }))
  }

  def filterRows(matrix: IndexedRowMatrix, samples: Array[String],
                 minMean: Double = Double.MinValue,
                 maxMean: Double = Double.MaxValue,
                 minSD: Double = 0,
                 maxSD: Double = Double.MaxValue): IndexedRowMatrix = {
    val broadcastSamples = SparkContext.getOrCreate().broadcast(samples)
    new IndexedRowMatrix(matrix.rows.filter(row => {
      val (mean, sd) = {
        val stats = meanAndVariance(MLibUtils.mllibVectorToDenseBreeze(row.vector))
        (stats.mean, Math.sqrt(stats.variance))
      }
      val toKeep = mean >= minMean && mean <= maxMean && sd >= minSD && sd <= maxSD
      if (!toKeep)
        // https://stackoverflow.com/questions/11940614/scala-and-slf4j-pass-multiple-parameters
        log.info("Excluding sample {} with mean={} and SD={}",
          Array[AnyRef](broadcastSamples.value(row.index.toInt), new java.lang.Double(mean), new java.lang.Double(sd)): _*)
      toKeep
    }))
  }

  def zscoreRows(matrix: IndexedRowMatrix): IndexedRowMatrix = ComputeZScores.time {
    new IndexedRowMatrix(matrix.rows.map(row => {
      val vector = MLibUtils.mllibVectorToDenseBreeze(row.vector)
      val stats = meanAndVariance(vector)
      IndexedRow(row.index, MLibUtils.breezeVectorToMLlib((vector - stats.mean) / Math.sqrt(stats.variance)))
    }))
  }

  def pcaNormalization(readMatrix: IndexedRowMatrix,
                       pveMeanFactor: Double = 0.7,
                       fixedToRemove: Option[Int] = None,
                       initialKFraction: Double = 0.10): IndexedRowMatrix = PCANormalization.time {
    readMatrix.rows.cache()

    // numRows returns maximum index, not actual number of rows
    val maxComponents = Math.min(readMatrix.rows.count(), readMatrix.numCols).toInt
    log.info("Maximum SVD components: {}", maxComponents)

    val (svd, bigK) = if (fixedToRemove.isDefined) {
      val bigK: Int = fixedToRemove.get
      log.info("Computing SVD with k={}", bigK)
      val svd = ComputeSVD.time {
        readMatrix.computeSVD(bigK, computeU = false)
      }
      (svd, bigK)
    } else {
      @tailrec def findK(k: Int): (SingularValueDecomposition[IndexedRowMatrix, org.apache.spark.mllib.linalg.Matrix], Int) = {
        // Compute SVD
        log.info("Computing SVD with k={}", k)
        val svd = ComputeSVD.time {
          readMatrix.computeSVD(k, computeU = false)
        }

        val componentVar = {
          val S = MLibUtils.mllibVectorToDenseBreeze(svd.s)
          S :* S
        }

        // estCutoff must be >= than actual cutoff, while minCutoff must be <= actual cutoff
        val minTotalVar = sum(componentVar)
        val estCutoff: Double = pveMeanFactor / maxComponents *
          (minTotalVar + Math.max(maxComponents - componentVar.length - 1, 0) * componentVar(-1))
        val minCutoff: Double = pveMeanFactor / maxComponents * minTotalVar

        @tailrec def findFirstLessThanCutoff(index: Int): Int = {
          if (index == componentVar.length || componentVar(index) < estCutoff)
            return index
          else
            return findFirstLessThanCutoff(index + 1)
        }
        val bigK = findFirstLessThanCutoff(0)

        // K is index of element one past the last component to be removed
        if (k == maxComponents || (bigK < componentVar.length && componentVar(bigK) < minCutoff)) {
          (svd, bigK)
        } else {
          log.info("SVD with k={} insufficient to determine K, recomputing.", k)
          findK(Math.min(2 * k, maxComponents))
        }
      }

      // Default mode where we determine K based on total variance
      findK(Math.min(Math.ceil(initialKFraction * maxComponents).toInt, maxComponents))
    }

    log.info("Removing top {} components in PCA normalization", bigK)
    val V = MLibUtils.mllibMatrixToDenseBreeze(svd.V)

    // Broadcast subset of V matrix
    val sc = SparkContext.getOrCreate()
    val C = sc.broadcast(V(::, 0 until bigK)) // Exclusive to K

    // Remove components with row centric approach
    new IndexedRowMatrix(readMatrix.rows.map(row => {
      val rd = MLibUtils.mllibVectorToDenseBreeze(row.vector)
      var rdStar = rd.copy
      for (col <- 0 until C.value.cols) {
        val component = C.value(::, col)
        rdStar :-= component * (component dot rd)
      }
      IndexedRow(row.index, MLibUtils.breezeVectorToMLlib(rdStar))
    }))

  }

  def normalizeReadDepth(readMatrix: ReadDepthMatrix,
                         minTargetMeanRD: Double = 10.0,
                         maxTargetMeanRD: Double = 500.0,
                         minSampleMeanRD: Double = 25.0,
                         maxSampleMeanRD: Double = 200.0,
                         minSampleSDRD: Double = 0.0,
                         maxSampleSDRD: Double = 150.0,
                         pveMeanFactor: Double = 0.7,
                         fixedToRemove: Option[Int] = None,
                         maxTargetSDRDStar: Double = 30.0,
                         initialKFraction: Double = 0.10): (IndexedRowMatrix, Array[ReferenceRegion]) = NormalizeReadDepths.time {

    // Filter I: Filter extreme targets and samples, then mean center the data
    val (centeredRdMatrix, targFilteredRdTargets) = ReadDepthFilterI.time {
      val (targFilteredrdMatrix, targFilteredRdTargets) = filterColumns(
        readMatrix.depth,
        readMatrix.targets,
        minMean = minTargetMeanRD,
        maxMean = maxTargetMeanRD)

      val sampFilteredRdMatrix = filterRows(
        targFilteredrdMatrix,
        readMatrix.samples,
        minMean = minSampleMeanRD,
        maxMean = maxSampleMeanRD,
        minSD = minSampleSDRD,
        maxSD = maxSampleSDRD)

      (meanCenterColumns(sampFilteredRdMatrix), targFilteredRdTargets)
    }

    // Can't use numRows to get row count as it is based on row indices
    centeredRdMatrix.rows.cache()
    log.info("After pre-normalization filter: {} samples by {} targets", centeredRdMatrix.rows.count, targFilteredRdTargets.length)

    // PCA normalization
    val rdStarMatrix = pcaNormalization(centeredRdMatrix,
      pveMeanFactor = pveMeanFactor,
      fixedToRemove = fixedToRemove,
      initialKFraction = initialKFraction)
    centeredRdMatrix.rows.unpersist() // Let Spark know centered data no longer needed

    // Filter II: Filter extremely variable targets
    val (targFilteredRdStarMatrix, targFilteredRdStarTargets) = ReadDepthFilterII.time {
      filterColumns(rdStarMatrix, targFilteredRdTargets, maxSD = maxTargetSDRDStar)
    }
    rdStarMatrix.rows.unpersist() // Let Spark know normalized data no longer needed

    // Z-score by sample
    val zMatrix = Normalization.zscoreRows(targFilteredRdStarMatrix)

    zMatrix.rows.cache()
    log.info("After post-normalization filter: {} samples by {} targets", zMatrix.rows.count, targFilteredRdStarTargets.length)
    (zMatrix, targFilteredRdStarTargets)
  }

}
