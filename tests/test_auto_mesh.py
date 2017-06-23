'''
Created on Oct 29, 2014

@author: opid30
'''

import os
import unittest
import tempfile
import json

import autoMesh

class Test(unittest.TestCase):

    def setUp(self):
        path = os.path.abspath(__file__)
        self.testDataDirectory = os.path.join(os.path.dirname(path), "data")


    def tes_checkForCorrelatedImages(self):
        testDataPath1 = os.path.join(self.testDataDirectory, "dictLoop_1.json")
        f = open(testDataPath1)
        dictLoop = json.loads(f.read())
        f.close()
        self.assertTrue(autoMesh.checkForCorrelatedImages(dictLoop))
        testDataPath2 = os.path.join(self.testDataDirectory, "dictLoop_2.json")
        f = open(testDataPath2)
        dictLoop = json.loads(f.read())
        f.close()
        self.assertFalse(autoMesh.checkForCorrelatedImages(dictLoop))



    def tes_autoMesh_identicalImages(self):
        beamline = "simulator"
        workingDir = tempfile.mkdtemp(prefix="autoMesh_", dir="/tmp_14_days/svensson")
        snapshotDir1 = "/scisoft/pxsoft/data/WORKFLOW_TEST_DATA/id30a1/snapshots/snapshots_sameimage_1"
        os.chmod(workingDir, 0755)
        (angleMinThickness, x1Pixels, y1Pixels, dxPixels, dyPixels, deltaPhiz, stdPhiz, imagePath, areTheSameImage) = \
            autoMesh.autoMesh(snapshotDir=snapshotDir1, workflowWorkingDir=workingDir, autoMeshWorkingDir=workingDir, \
                          loopMaxWidth=330, loopMinWidth=250, prefix="snapshot", debug=False, findLargestMesh=False)
        self.assertTrue(areTheSameImage)
        snapshotDir2 = "/scisoft/pxsoft/data/WORKFLOW_TEST_DATA/id30a1/snapshots/snapshots_20151013-152805_jhp9cH"
        (angleMinThickness, x1Pixels, y1Pixels, dxPixels, dyPixels, deltaPhiz, stdPhiz, imagePath, areTheSameImage) = \
            autoMesh.autoMesh(snapshotDir=snapshotDir2, workflowWorkingDir=workingDir, autoMeshWorkingDir=workingDir, \
                          loopMaxWidth=330, loopMinWidth=250, prefix="snapshot", debug=False, findLargestMesh=False)
        self.assertFalse(areTheSameImage)
        snapshotDir3 = "/scisoft/pxsoft/data/WORKFLOW_TEST_DATA/id30a1/snapshots/snapshots_20151211-133113_unt8EX"
        (angleMinThickness, x1Pixels, y1Pixels, dxPixels, dyPixels, deltaPhiz, stdPhiz, imagePath, areTheSameImage) = \
            autoMesh.autoMesh(snapshotDir=snapshotDir3, workflowWorkingDir=workingDir, autoMeshWorkingDir=workingDir, \
                          loopMaxWidth=330, loopMinWidth=250, prefix="snapshot", debug=False, findLargestMesh=False)
        self.assertFalse(areTheSameImage)
        snapshotDir4 = "/scisoft/pxsoft/data/WORKFLOW_TEST_DATA/id30a1/snapshots//snapshots_20160203-092805_rJIRxO"
        (angleMinThickness, x1Pixels, y1Pixels, dxPixels, dyPixels, deltaPhiz, stdPhiz, imagePath, areTheSameImage) = \
            autoMesh.autoMesh(snapshotDir=snapshotDir4, workflowWorkingDir=workingDir, autoMeshWorkingDir=workingDir, \
                          loopMaxWidth=330, loopMinWidth=250, prefix="snapshot", debug=False, findLargestMesh=False)
        self.assertTrue(areTheSameImage)

    def test_autoMesh(self):
        beamline = "simulator"
        snapshotDir = os.path.join(self.testDataDirectory, "snapshots_20141128-084026")
        workingDir = tempfile.mkdtemp(prefix="autoMesh_", dir="/tmp_14_days/svensson")
        print(workingDir)
        os.chmod(workingDir, 0755)
        identifier = snapshotDir[-22:-7]
        (newPhi, x1Pixels, y1Pixels, dxPixels, dyPixels, deltaPhizPixels, stdPhizPixels, imagePath, areTheSameImage) = \
            autoMesh.autoMesh(snapshotDir, workingDir, workingDir,
                    debug=True,
                    loopMaxWidth=0.8 * 420,
                    loopMinWidth=0.35 * 420,
                    findLargestMesh=False)
        print("Grid coordinates in pixels: x1=%.1f y1=%.1f dx=%.1f dy=%.1f" % (x1Pixels, y1Pixels, dxPixels, dyPixels))
        PIXELS_PER_MM = 420.0
        deltaPhizPixels = 0.0
        xOffsetLeft = 0.05
        xOffsetRight = 0.05
        beamSize = 0.1
        overSamplingX = 1.5
        overSamplingY = 1.5
        x1 = x1Pixels / PIXELS_PER_MM - xOffsetLeft
        y1 = -(y1Pixels + dyPixels + deltaPhizPixels) / PIXELS_PER_MM
        dx_mm = dxPixels / PIXELS_PER_MM + xOffsetLeft + xOffsetRight
        dy_mm = dyPixels / PIXELS_PER_MM
        steps_x = int((dxPixels / PIXELS_PER_MM + xOffsetLeft + xOffsetRight) / beamSize * overSamplingX)
        steps_y = int(dyPixels / PIXELS_PER_MM / beamSize * overSamplingY)
        # Check that we have at least three lines:
        if steps_y == 1 or steps_y == 2:
            steps_y = 3
            y1 -= dyPixels / PIXELS_PER_MM / 4.0
            dy_mm *= 1.5
        grid_info = {"x1": x1,
                     "y1": y1,
                     "dx_mm": dx_mm,
                     "dy_mm": dy_mm,
                     "steps_x": steps_x,
                     "steps_y": steps_y}
        print("Auto grid_info: %r" % grid_info)
        print workingDir
        resultImagePath = os.path.join(workingDir, "snapshot_automesh.png")
        autoMesh.plotMesh(beamline, imagePath, grid_info, PIXELS_PER_MM, workingDir)
        os.system("display %s" % resultImagePath)
