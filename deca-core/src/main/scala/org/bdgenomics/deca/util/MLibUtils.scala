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
package org.bdgenomics.deca.util

import breeze.linalg.{ DenseMatrix, DenseVector, SliceVector, SparseVector }

/**
 * Adapted from https://github.com/amplab/keystone/blob/master/src/main/scala/utils/MLlibUtils.scala
 */
object MLibUtils {

  /** Convert an MLlib vector to a Breeze dense vector */
  def mllibVectorToDenseBreeze(vector: org.apache.spark.mllib.linalg.Vector): DenseVector[Double] = {
    vector match {
      case dense: org.apache.spark.mllib.linalg.DenseVector => new DenseVector[Double](dense.values)
      case _ => new DenseVector[Double](vector.toArray)
    }
  }

  /** Convert an MLlib matrix to a Breeze dense matrix */
  def mllibMatrixToDenseBreeze(matrix: org.apache.spark.mllib.linalg.distributed.IndexedRowMatrix): DenseMatrix[Double] = {
    // Breeze is column ordered so create transposed matrix from rows
    val breezeMatrix = new DenseMatrix[Double](matrix.numCols.toInt, matrix.numRows.toInt)
    matrix.rows.collect.foreach(row => {
      breezeMatrix(::, row.index.toInt) := mllibVectorToDenseBreeze(row.vector)
    })
    breezeMatrix.t
  }

  def mllibMatrixToDenseBreeze(matrix: org.apache.spark.mllib.linalg.Matrix): DenseMatrix[Double] = {
    matrix match {
      case dense: org.apache.spark.mllib.linalg.DenseMatrix => {
        if (!dense.isTransposed) {
          new DenseMatrix[Double](dense.numRows, dense.numCols, dense.values)
        } else {
          val breezeMatrix = new DenseMatrix[Double](dense.numRows, dense.numCols, dense.values)
          breezeMatrix.t
        }
      }

      case _ => new DenseMatrix[Double](matrix.numRows, matrix.numCols, matrix.toArray)
    }
  }

  /** Convert a Breeze vector to an MLlib vector, maintaining underlying data structure (sparse vs dense) */
  def breezeVectorToMLlib(breezeVector: breeze.linalg.Vector[Double]): org.apache.spark.mllib.linalg.Vector = {
    breezeVector match {
      case v: SliceVector[Int, Double] => breezeVectorToMLlib(v.toDenseVector)
      case v: DenseVector[Double] =>
        if (v.offset == 0 && v.stride == 1 && v.length == v.data.length) {
          new org.apache.spark.mllib.linalg.DenseVector(v.data)
        } else {
          new org.apache.spark.mllib.linalg.DenseVector(v.toArray) // Can't use underlying array directly, so make a new one
        }
      case v: SparseVector[Double] =>
        if (v.index.length == v.used) {
          new org.apache.spark.mllib.linalg.SparseVector(v.length, v.index, v.data)
        } else {
          new org.apache.spark.mllib.linalg.SparseVector(v.length, v.index.slice(0, v.used), v.data.slice(0, v.used))
        }
      case v: breeze.linalg.Vector[_] =>
        sys.error("Unsupported Breeze vector type: " + v.getClass.getName)
    }
  }
}
