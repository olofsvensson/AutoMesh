# coding: utf-8
# /*##########################################################################
# Copyright (C) 2017 European Synchrotron Radiation Facility
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ############################################################################*/

import os
import shutil
import pathlib
import unittest
import tempfile

import lib_auto_mesh



class Test(unittest.TestCase):
    def setUp(self):
        self.test_data_directory = pathlib.Path(__file__).parent / "data"
        self.working_dir = tempfile.mkdtemp(prefix="autoMesh_")
        os.chmod(self.working_dir, 0o755)

    def tearDown(self) -> None:
        shutil.rmtree(self.working_dir)

    def test_tungsten(self):
        snapshot_dir = self.test_data_directory / "tungsten"
        print(snapshot_dir)
        pixelPerMm = 1024
        print("PixelPerMm: {0}".format(pixelPerMm))
        tungstenPixels = 0.020 * pixelPerMm / 2
        print("Tungsten pixels: {0}".format(tungstenPixels))
        prefix = "snapshot"
        background_image = os.path.join(snapshot_dir, "%s_background.png" % prefix)
        print(self.working_dir)
        deltaX, deltaY, deltaZ = lib_auto_mesh.findDeltaToCentre(
            snapshot_dir,
            self.working_dir,
            loop_width=tungstenPixels,
            debug=True
        )
        print(deltaX, deltaY, deltaZ)
        print(deltaX / pixelPerMm, deltaY / pixelPerMm)


if __name__ == "__main__":
    unittest.main()
