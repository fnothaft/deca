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

import java.io.{ InputStream, OutputStream }

import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.fs.{ FileSystem, FileUtil, Path }
import org.apache.hadoop.io.IOUtils
import org.bdgenomics.utils.misc.Logging

/**
 * Helper object to merge sharded files together adapted from ADAM
 */
object FileMerger extends Logging {

  /**
   * Merges together sharded files.
   *
   * @param fs The file system implementation to use.
   * @param outputPath The location to write the merged file at.
   * @param bodyPath The location where the sharded files have been written.
   * @param optHeaderPath Optionally, the location where a header file has
   *   been written.
   */
  private[deca] def mergeFiles(conf: Configuration,
                               fs: FileSystem,
                               outputPath: Path,
                               bodyPath: Path,
                               optHeaderPath: Option[Path] = None) {

    val os = fs.create(outputPath)

    def copy(inputPath: Path, los: OutputStream) {
      val lis = fs.open(inputPath)
      try {
        IOUtils.copyBytes(lis, los, conf, false)
      } finally {
        lis.close()
      }
    }

    // optionally copy the header
    optHeaderPath.foreach(headerPath => {
      log.info("Copying header file (%s)".format(headerPath))
      copy(headerPath, os)
    })

    // get a list of all of the files in the tail file
    val bodyFiles = fs.globStatus(new Path("%s/part-*".format(bodyPath)))
      .toSeq
      .map(_.getPath)
      .sortBy(_.getName)

    val numFiles = bodyFiles.length
    var filesCopied = 1
    bodyFiles.foreach(file => {
      log.info("Copying file %s, file %d of %d.".format(file.toString, filesCopied, numFiles))
      copy(file, os)
      filesCopied += 1
    })

    os.flush()
    os.close()

    // Cleanup temporary files
    optHeaderPath.foreach(headPath => fs.delete(headPath, true))
    fs.delete(bodyPath, true)
  }

}
