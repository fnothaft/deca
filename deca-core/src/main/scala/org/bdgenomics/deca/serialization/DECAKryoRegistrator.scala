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
package org.bdgenomics.deca.serialization

import com.esotericsoftware.kryo.Kryo
import org.apache.spark.serializer.KryoRegistrator
import org.bdgenomics.adam.serialization.ADAMKryoRegistrator

/**
 * Created by mlinderman on 6/16/17.
 */

class DECAKryoRegistrator extends KryoRegistrator {

  private val akr = new ADAMKryoRegistrator()

  override def registerClasses(kryo: Kryo) {

    // Register ADAM's requirements
    akr.registerClasses(kryo)

    // breeze.linalg
    kryo.register(classOf[breeze.linalg.DenseVector[Double]])
    kryo.register(classOf[breeze.linalg.DenseVector[Boolean]])
    kryo.register(classOf[breeze.linalg.DenseVector$mcD$sp])
    kryo.register(classOf[breeze.linalg.BitVector])
    kryo.register(classOf[breeze.linalg.DenseMatrix$mcD$sp])

    // org.apache.spark.mllib
    kryo.register(classOf[org.apache.spark.mllib.linalg.DenseMatrix])
    kryo.register(classOf[org.apache.spark.mllib.linalg.DenseVector])
    kryo.register(classOf[org.apache.spark.mllib.linalg.distributed.IndexedRow])
    kryo.register(classOf[scala.Array[org.apache.spark.mllib.linalg.distributed.IndexedRow]])
    kryo.register(classOf[org.apache.spark.mllib.stat.MultivariateOnlineSummarizer])

    // org.bdgenomics.deca
    kryo.register(classOf[org.bdgenomics.deca.hmm.TransitionProbabilities])

    // scala
    kryo.register(classOf[scala.Array[Boolean]])
    kryo.register(classOf[scala.Array[Double]])
  }
}
