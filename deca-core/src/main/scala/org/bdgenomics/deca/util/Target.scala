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

import org.bdgenomics.adam.models.ReferenceRegion

/**
 * Created by mlinderman on 6/5/17.
 */
object Target {

  /**
   * Generates ReferenceRegion from compact, inclusive region description
   *
   * @param region String of "contig:start-end"
   * @return Equivalent ReferenceRegion
   */
  def regionToReferenceRegion(region: String): ReferenceRegion = {
    val fields = region.split(Array(':', '-'))
    new ReferenceRegion(fields(0), fields(1).toLong - 1, fields(2).toLong)
  }

}
