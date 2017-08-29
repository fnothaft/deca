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

import breeze.numerics.abs
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.deca.util.MLibUtils

/**
 * Created by mlinderman on 3/9/17.
 */
class CoverageSuite extends DecaFunSuite {

  sparkTest("Finds correct first target") {
    val inputTargets = resourceUrl("EXOME.interval_list")
    val features = sc.loadFeatures(inputTargets.toString)

    val targetsDriver = features.rdd
      .map(ReferenceRegion.unstranded(_))
      .collect
      .zipWithIndex
      .sortBy(_._1)

    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 16448862, 16448952)) === Some(0))
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 16448800, 16448810)) === Some(0))
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 20074000, 20074010)).isEmpty)
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 17071000, 17071010)) === Some(4))
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 16449200, 16449300)) === Some(1))
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 16448824, 16449024)) === Some(0))
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 16448834, 16449024)) === Some(0))
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 16448804, 16449000)) === Some(0))
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("21", 16448804, 16449000)) === Some(0))
    assert(Coverage.findFirstTarget(targetsDriver, ReferenceRegion("22", 19127400, 19127500)) === Some(156))
  }

  sparkTest("computes coverage for targets") {
    val inputBam = resourceUrl("HG00107_target1.bam")
    val reads = sc.loadAlignments(inputBam.toString)

    val inputTargets = resourceUrl("EXOME.interval_list")
    val features = sc.loadFeatures(inputTargets.toString)

    val depths = Coverage.coverageMatrix(Seq(reads), features, minMapQ = 20)
    assert(depths.numSamples() === 1 && depths.numTargets() === 300)

    val matrix = MLibUtils.mllibMatrixToDenseBreeze(depths.depth)

    assert(abs(matrix(0, 0) - 18.31) < 0.005) // Should be target 22:16448824-16449023
    assert(abs(matrix(0, 1) - 58.86) < 0.005) // There are reads overlapping next target
    assert(matrix(0, 2) === 0.0)
  }

  sparkTest("counts fragments not just reads") {
    val inputBam = resourceUrl("HG00146_target51.bam")
    val reads = sc.loadAlignments(inputBam.toString)

    val inputTargets = resourceUrl("EXOME.interval_list")
    val features = sc.loadFeatures(inputTargets.toString)

    val depths = Coverage.coverageMatrix(Seq(reads), features, minMapQ = 20)
    assert(depths.numSamples() === 1 && depths.numTargets() === 300)

    val matrix = MLibUtils.mllibMatrixToDenseBreeze(depths.depth)

    assert(abs(matrix(0, 50) - 407.32) < 0.005) // Should be target 22:17662374-17662466
  }

  sparkTest("counts imperfectly paired fragments") {
    val inputBam = resourceUrl("HG00116_column.bam")
    val reads = sc.loadAlignments(inputBam.toString)

    val inputTargets = resourceUrl("HG00116_column.interval_list")
    val features = sc.loadFeatures(inputTargets.toString)

    val depths = Coverage.coverageMatrix(Seq(reads), features, minMapQ = 20)
    assert(depths.numSamples() === 1 && depths.numTargets() === 1)

    val matrix = MLibUtils.mllibMatrixToDenseBreeze(depths.depth)

    assert(abs(matrix(0, 0) - 343) < 1e-8) // Should be target 22:19398239
  }

  sparkTest("incorporates CIGAR into fragment analysis") {
    val inputBam = resourceUrl("HG00146_overlap.bam")
    val reads = sc.loadAlignments(inputBam.toString)

    val inputTargets = resourceUrl("EXOME.interval_list")
    val features = sc.loadFeatures(inputTargets.toString)

    val depths = Coverage.coverageMatrix(Seq(reads), features, minMapQ = 20)
    assert(depths.numSamples() === 1 && depths.numTargets() === 300)

    val matrix = MLibUtils.mllibMatrixToDenseBreeze(depths.depth)

    assert(abs(matrix(0, 34) - 0.31) < 0.005) // 22:17600882-17601080
    assert(abs(matrix(0, 35) - 1.00) < 0.005) // 22:17601082-17601091

  }
}
