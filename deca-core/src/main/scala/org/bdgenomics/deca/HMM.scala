package org.bdgenomics.deca

import org.apache.spark.SparkContext
import org.apache.spark.storage.StorageLevel
import org.bdgenomics.adam.rdd.feature.FeatureRDD
import org.bdgenomics.deca.Timers._
import org.bdgenomics.deca.coverage.ReadDepthMatrix
import org.bdgenomics.deca.hmm.{ SampleModel, TransitionProbabilities }
import org.bdgenomics.deca.util.MLibUtils
import org.bdgenomics.formats.avro.{ Feature, Strand }
import org.bdgenomics.utils.misc.Logging

/**
 * Created by mlinderman on 4/4/17.
 */
object HMM extends Serializable with Logging {

  def discoverCNVs(readMatrix: ReadDepthMatrix,
                   M: Double = 3, T: Double = 6, p: Double = 1e-8, D: Double = 70000,
                   minSomeQuality: Double = 30.0): FeatureRDD = DiscoverCNVs.time {

    val sc = SparkContext.getOrCreate()

    // Generate transition probabilities from targets
    val transProb = sc.broadcast(TransitionProbabilities(readMatrix.targets, D = D, p = p, q = 1.0 / T))

    // Broadcast samples and targets so that features can be fixed up
    val targets = sc.broadcast(readMatrix.targets)
    val samples = sc.broadcast(readMatrix.samples)

    val cnvs = readMatrix.depth.rows.flatMap(obs => {
      // Create per-sample HMM model
      val model = SampleModel(obs.vector, transProb.value, M, p)

      // Discover CNVs
      val per_sample_cnvs = model.discoverCNVs(minSomeQuality)

      // Refine feature descriptions with coordinates, sample, etc.
      per_sample_cnvs.map(raw_feature => {
        val attr = raw_feature.getAttributes

        val start_index = attr.get("START_TARGET").toInt
        val end_index = attr.get("END_TARGET").toInt

        val start_target = targets.value(start_index)
        val end_target = targets.value(end_index)

        // TODO: Filter out contig spanning CNVs

        val builder = Feature.newBuilder(raw_feature)
        builder.setSource(samples.value(obs.index.toInt))
        builder.setContigName(start_target.referenceName)
        builder.setStart(start_target.start)
        builder.setEnd(end_target.end)
        builder.setStrand(Strand.INDEPENDENT)

        // Transform START_TARGET and END_TARGET to be 1-indexed
        attr.put("START_TARGET", (start_index + 1).toString)
        attr.put("END_TARGET", (end_index + 1).toString)

        builder.setAttributes(attr)

        builder.build()
      })
    })

    FeatureRDD.inferSequenceDictionary(cnvs, optStorageLevel = Some(StorageLevel.MEMORY_ONLY_SER))
  }
}
