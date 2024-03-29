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
import json
import shutil
import unittest
import tempfile

import lib_auto_mesh

SCISOFT_DIR = "/scisoft/pxsoft/data/WORKFLOW_TEST_DATA/id30a1/snapshots"


class Test(unittest.TestCase):
    def setUp(self):
        path = os.path.abspath(__file__)
        self.test_data_directory = os.path.join(os.path.dirname(path), "data")
        self.working_dir = tempfile.mkdtemp(prefix="autoMesh_")
        os.chmod(self.working_dir, 0o755)

    def tearDown(self) -> None:
        shutil.rmtree(self.working_dir)

    def test_checkForCorrelatedImages(self):
        test_data_path1 = os.path.join(self.test_data_directory, "dictLoop_1.json")
        f = open(test_data_path1)
        dict_loop = json.loads(f.read())
        f.close()
        self.assertTrue(lib_auto_mesh.checkForCorrelatedImages(dict_loop))
        test_data_path2 = os.path.join(self.test_data_directory, "dictLoop_2.json")
        f = open(test_data_path2)
        dict_loop = json.loads(f.read())
        f.close()
        self.assertFalse(lib_auto_mesh.checkForCorrelatedImages(dict_loop))

    def test_autoMesh_identicalImages_1(self):
        snapshot_dir1 = os.path.join(self.test_data_directory, "snapshots_sameimage_1")
        (
            angle_min_thickness,
            x1_pixels,
            y1_pixels,
            dx_pixels,
            dy_pixels,
            delta_phiz,
            std_phiz,
            image_path,
            are_the_same_image,
        ) = lib_auto_mesh.autoMesh(
            snapshot_dir=snapshot_dir1,
            workflow_working_dir=self.working_dir,
            auto_mesh_working_dir=self.working_dir,
            loop_max_width=330,
            loop_min_width=250,
            prefix="snapshot",
            debug=False,
            find_largest_mesh=False,
        )
        self.assertTrue(are_the_same_image)

    @unittest.skipIf(
        not os.path.exists(
            os.path.join(SCISOFT_DIR, "snapshots_20151013-152805_jhp9cH")
        ),
        "Test directory missing",
    )
    def test_autoMesh_identicalImages_2(self):
        snapshot_dir2 = os.path.join(SCISOFT_DIR, "snapshots_20151013-152805_jhp9cH")
        (
            angle_min_thickness,
            x1_pixels,
            y1_pixels,
            dx_pixels,
            dy_pixels,
            delta_phiz,
            std_phiz,
            image_path,
            are_the_same_image,
        ) = lib_auto_mesh.autoMesh(
            snapshot_dir=snapshot_dir2,
            workflow_working_dir=self.working_dir,
            auto_mesh_working_dir=self.working_dir,
            loop_max_width=330,
            loop_min_width=250,
            prefix="snapshot",
            debug=False,
            find_largest_mesh=False,
        )
        self.assertFalse(are_the_same_image)

    @unittest.skipIf(
        not os.path.exists(
            os.path.join(SCISOFT_DIR, "snapshots_20151211-133113_unt8EX")
        ),
        "Test directory missing",
    )
    def test_autoMesh_identicalImages_3(self):
        snapshot_dir3 = os.path.join(SCISOFT_DIR, "snapshots_20151211-133113_unt8EX")
        (
            angle_min_thickness,
            x1_pixels,
            y1_pixels,
            dx_pixels,
            dy_pixels,
            delta_phiz,
            std_phiz,
            image_path,
            are_the_same_image,
        ) = lib_auto_mesh.autoMesh(
            snapshot_dir=snapshot_dir3,
            workflow_working_dir=self.working_dir,
            auto_mesh_working_dir=self.working_dir,
            loop_max_width=330,
            loop_min_width=250,
            prefix="snapshot",
            debug=False,
            find_largest_mesh=False,
        )
        self.assertFalse(are_the_same_image)

    @unittest.skipIf(
        not os.path.exists(
            os.path.join(SCISOFT_DIR, "snapshots_20160203-092805_rJIRxO")
        ),
        "Test directory missing",
    )
    def test_autoMesh_identicalImages_4(self):
        snapshot_dir4 = os.path.join(SCISOFT_DIR, "snapshots_20160203-092805_rJIRxO")
        (
            angle_min_thickness,
            x1_pixels,
            y1_pixels,
            dx_pixels,
            dy_pixels,
            delta_phiz,
            std_phiz,
            image_path,
            are_the_same_image,
        ) = lib_auto_mesh.autoMesh(
            snapshot_dir=snapshot_dir4,
            workflow_working_dir=self.working_dir,
            auto_mesh_working_dir=self.working_dir,
            loop_max_width=330,
            loop_min_width=250,
            prefix="snapshot",
            debug=False,
            find_largest_mesh=False,
        )
        self.assertTrue(are_the_same_image)

    def test_autoMesh(self):
        snapshot_dir = os.path.join(
            self.test_data_directory, "snapshots_20141128-084026"
        )
        (
            new_phi,
            x1_pixels,
            y1_pixels,
            dx_pixels,
            dy_pixels,
            delta_phiz_pixels,
            std_phiz_pixels,
            image_path,
            are_the_same_image,
        ) = lib_auto_mesh.autoMesh(
            snapshot_dir,
            self.working_dir,
            self.working_dir,
            debug=False,
            loop_max_width=0.35 * 608,
            loop_min_width=0.5 * 608,
            find_largest_mesh=True,
        )
        print(
            "Grid coordinates in pixels: x1=%.1f y1=%.1f dx=%.1f dy=%.1f"
            % (x1_pixels, y1_pixels, dx_pixels, dy_pixels)
        )
        PIXELS_PER_MM = 420.0
        delta_phiz_pixels = 0.0
        x_offset_left = 0.05
        x_offset_right = 0.05
        beam_size = 0.1
        over_sampling_x = 1.5
        over_sampling_y = 1.5
        x1 = x1_pixels / PIXELS_PER_MM - x_offset_left
        y1 = -(y1_pixels + dy_pixels + delta_phiz_pixels) / PIXELS_PER_MM
        dx_mm = dx_pixels / PIXELS_PER_MM + x_offset_left + x_offset_right
        dy_mm = dy_pixels / PIXELS_PER_MM
        steps_x = int(
            (dx_pixels / PIXELS_PER_MM + x_offset_left + x_offset_right)
            / beam_size
            * over_sampling_x
        )
        steps_y = int(dy_pixels / PIXELS_PER_MM / beam_size * over_sampling_y)
        # Check that we have at least three lines:
        if steps_y == 1 or steps_y == 2:
            steps_y = 3
            y1 -= dy_pixels / PIXELS_PER_MM / 4.0
            dy_mm *= 1.5
        grid_info = {
            "x1": x1,
            "y1": y1,
            "dx_mm": dx_mm,
            "dy_mm": dy_mm,
            "steps_x": steps_x,
            "steps_y": steps_y,
        }
        print("Auto grid_info: %r" % grid_info)
        result_image_path = os.path.join(self.working_dir, "snapshot_automesh.png")
        lib_auto_mesh.plotMesh(
            image_path, grid_info, PIXELS_PER_MM, self.working_dir, show_plot=False
        )
        # os.system("display %s" % result_image_path)
        self.assertTrue(os.path.exists(result_image_path))


if __name__ == "__main__":
    unittest.main()
