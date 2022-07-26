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
"""
Workflow library module for grid / mesh
"""

__author__ = "Olof Svensson"
__contact__ = "svensson@esrf.eu"
__copyright__ = "ESRF, 2017"
__updated__ = "2017-06-23"


import os
import numpy
import logging
import imageio
import scipy.ndimage

# import matplotlib
# matplotlib.use("Agg", warn=False)
import matplotlib.pyplot as pyplot
import pylab

logging.basicConfig(level=logging.INFO)


def autoMesh(
    snapshot_dir,
    workflow_working_dir,
    auto_mesh_working_dir,
    loop_max_width=300,
    loop_min_width=150,
    prefix="snapshot",
    debug=False,
    find_largest_mesh=False,
):
    angle_min_thickness = None
    x1_pixels = None
    y1_pixels = None
    dx_pixels = None
    dy_pixels = None
    delta_phiz = None
    std_phiz = None
    image_path = None
    os.chmod(auto_mesh_working_dir, 0o755)
    background_image = os.path.join(snapshot_dir, "%s_background.png" % prefix)
    dict_loop = {}
    for omega in [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]:
        logging.info("Analysing snapshot image at omega = %d degrees" % omega)
        image_path = os.path.join(snapshot_dir, "%s_%03d.png" % (prefix, omega))
        raw_img = readImage(image_path)

        if debug:
            plot_img(
                raw_img,
                os.path.join(auto_mesh_working_dir, "rawImage_%03d.png" % omega),
            )

        if background_image.endswith(".npy"):
            background = numpy.load(background_image)
        else:
            background = imageio.imread(background_image, as_gray=True)

        difference_image = numpy.abs(background - raw_img)

        if debug:
            plot_img(
                difference_image,
                plot_path=os.path.join(
                    auto_mesh_working_dir, "differenceImage_%03d.png" % omega
                ),
            )

        #        print numpyImage
        filteredImage = filterDifferenceImage(difference_image)
        if debug:
            plot_img(
                filteredImage,
                plot_path=os.path.join(
                    auto_mesh_working_dir, "filteredImage_%03d.png" % omega
                ),
            )
        (list_index, list_upper, list_lower) = loopExam(filteredImage)
        if debug:
            pylab.plot(list_index, list_upper, "+")
            pylab.plot(list_index, list_lower, "+")
            raw_img = imageio.imread(image_path, as_gray=True)
            imgshape = raw_img.shape
            extent = (0, imgshape[1], 0, imgshape[0])
            pylab.axes(extent)
            pylab.savefig(
                os.path.join(auto_mesh_working_dir, "shapePlot_%03d.png" % omega)
            )
            pyplot.close()
        dict_loop["%d" % omega] = (list_index, list_upper, list_lower)
    are_the_same_image = checkForCorrelatedImages(dict_loop)
    if not are_the_same_image:
        image000_path = os.path.join(snapshot_dir, "%s_%03d.png" % (prefix, 0))
        ny, nx = imageio.imread(image000_path, as_gray=True).shape
        (
            angle_min_thickness,
            x1_pixels,
            y1_pixels,
            dx_pixels,
            dy_pixels,
            delta_phiz,
            std_phiz,
        ) = findOptimalMesh(
            dict_loop,
            snapshot_dir,
            nx,
            ny,
            auto_mesh_working_dir,
            debug=debug,
            loop_max_width=loop_max_width,
            loop_min_width=loop_min_width,
            find_largest_mesh=find_largest_mesh,
        )
        image_path = None
        if angle_min_thickness is not None:
            image_path = os.path.join(
                snapshot_dir, "%s_%03d.png" % (prefix, angle_min_thickness)
            )
    return (
        angle_min_thickness,
        x1_pixels,
        y1_pixels,
        dx_pixels,
        dy_pixels,
        delta_phiz,
        std_phiz,
        image_path,
        are_the_same_image,
    )


def subtractBackground(image, background_image):
    return numpy.abs(image - background_image)


def filterDifferenceImage(difference_image, threshold_value=30):
    """
    First applies a threshold of default value 30.
    Then erodes the image twice, and then dilates the image twice.
    """
    threshold_value = 30
    binary_image = difference_image >= threshold_value
    filtered_image = scipy.ndimage.morphology.binary_erosion(binary_image).astype(
        binary_image.dtype
    )
    filtered_image = scipy.ndimage.morphology.binary_erosion(filtered_image).astype(
        binary_image.dtype
    )
    filtered_image = scipy.ndimage.morphology.binary_dilation(filtered_image).astype(
        binary_image.dtype
    )
    filtered_image = scipy.ndimage.morphology.binary_dilation(filtered_image).astype(
        binary_image.dtype
    )
    return filtered_image


def loopExam(filtered_image):
    """
    This method examines the loop in one image.

    """
    ny, nx = filtered_image.shape
    shape_list_index = []
    shape_list_upper = []
    shape_list_lower = []
    for index_x in range(nx):
        column = filtered_image[:, index_x]
        indices = numpy.where(column)[0]
        if len(indices) > 0:
            shape_list_index.append(index_x)
            shape_list_upper.append(indices[0])
            shape_list_lower.append(indices[-1])
    array_index = numpy.array(shape_list_index)
    array_upper = ny - numpy.array(shape_list_upper)
    array_lower = ny - numpy.array(shape_list_lower)
    return (array_index.tolist(), array_upper.tolist(), array_lower.tolist())


def checkForCorrelatedImages(dict_loop):
    # Check if all the indices are the same
    first_list_index = None
    first_upper_index = None
    first_lower_index = None
    are_the_same = True
    for loop_index in dict_loop:
        list_index, list_upper, list_lower = dict_loop[loop_index]
        if first_list_index is None:
            first_list_index = list_index
            first_upper_index = list_upper
            first_lower_index = list_lower
        elif cmp(first_list_index, list_index):
            are_the_same = False
            break
        elif cmp(first_upper_index, list_upper):
            are_the_same = False
            break
        elif cmp(first_lower_index, list_lower):
            are_the_same = False
            break
    return are_the_same


def findOptimalMesh(
    dict_loop,
    snapshot_dir,
    nx,
    ny,
    auto_mesh_working_dir,
    loop_max_width=300,
    loop_min_width=150,
    debug=False,
    find_largest_mesh=False,
):
    array_phiz = None
    min_thickness = None
    angle_min_thickness = None
    max_thickness = None
    angle_max_thickness = None
    mesh_xmax = None
    std_phiz = None
    loop_min_width = int(loop_min_width)
    loop_max_width = int(loop_max_width)
    for omega in [0, 30, 60, 90, 120, 150]:
        str_omega1 = "%d" % omega
        str_omega2 = "%03d" % (omega + 180)
        (list_index1, list_upper1, list_lower1) = dict_loop[str_omega1]
        (list_index2, list_upper2, list_lower2) = dict_loop[str_omega2]
        array_index1 = numpy.array(list_index1)
        array_index2 = numpy.array(list_index2)
        if len(array_index1) == 0 or len(array_index2) == 0:
            break
        if mesh_xmax is None:
            mesh_xmax = numpy.max(array_index1)
        elif mesh_xmax > numpy.max(array_index1):
            mesh_xmax = numpy.max(array_index1)
        if mesh_xmax > numpy.max(array_index2):
            mesh_xmax = numpy.max(array_index2)
    if mesh_xmax is None:
        mesh_xmax = loop_max_width
    mesh_xmin = mesh_xmax - loop_max_width
    print(mesh_xmin, mesh_xmax)
    n_phi = 0
    for omega in [0, 30, 60, 90, 120, 150]:
        str_omega1 = "%d" % omega
        str_omega2 = "%03d" % (omega + 180)
        (array_index1, array_upper1, array_lower1) = dict_loop[str_omega1]
        (array_index2, array_upper2, array_lower2) = dict_loop[str_omega2]
        array_index = numpy.arange(mesh_xmin, mesh_xmax - 1)
        # Exclude region in horizontal center
        x_exclud_min = numpy.float64(nx / 2 - 20)
        x_exclud_max = numpy.float64(nx / 2 + 20)
        indices1 = numpy.where(
            ((array_index1 > mesh_xmin) & (array_index1 < x_exclud_min))
            | ((array_index1 < mesh_xmax) & (array_index1 > x_exclud_max))
        )
        indices2 = numpy.where(
            ((array_index2 > mesh_xmin) & (array_index2 < x_exclud_min))
            | ((array_index2 < mesh_xmax) & (array_index2 > x_exclud_max))
        )
        array_upper1 = numpy.array(array_upper1)[indices1]
        array_upper2 = numpy.array(array_upper2)[indices2]
        array_lower1 = numpy.array(array_lower1)[indices1]
        array_lower2 = numpy.array(array_lower2)[indices2]
        if True:
            pylab.plot(array_upper1, "+", color="red")
            pylab.plot(array_upper2, "+", color="blue")
            pylab.plot(array_lower1, "+", color="red")
            pylab.plot(array_lower2, "+", color="blue")
            pylab.savefig(
                os.path.join(
                    auto_mesh_working_dir,
                    "shapePlot_%03d_%03d.png" % (omega, omega + 180),
                )
            )
            pylab.close()
        phiz = None
        if (array_upper1.shape == array_lower2.shape) and (
            array_upper2.shape == array_lower1.shape
        ):
            phiz1 = (array_upper1 + array_lower2) / 2.0
            phiz2 = (array_upper2 + array_lower1) / 2.0
            phiz = (phiz1 + phiz2) / 2.0
        elif array_upper1.shape == array_lower2.shape:
            phiz = (array_upper1 + array_lower2) / 2.0
        elif array_upper2.shape == array_lower1.shape:
            phiz = (array_upper2 + array_lower1) / 2.0
        if phiz is not None:
            if debug:
                pylab.plot(phiz, "+")
                pylab.savefig(
                    os.path.join(
                        auto_mesh_working_dir,
                        "phiz_%03d_%03d.png" % (omega, omega + 180),
                    )
                )
                pylab.close()
            if array_phiz is None:
                array_phiz = numpy.array(phiz)
                n_phi += 1
            else:
                if array_phiz.shape == phiz.shape:
                    array_phiz += phiz
                    n_phi += 1
    if n_phi > 0:
        array_phiz /= n_phi
        # Cut off 20 points from each side of array
        # in order to remove artifacts at end points
        array_phiz = array_phiz[20:-20]
        # arrayIndexPhiz = array_index[20:-20]
        if True:
            pylab.plot(array_phiz, "+")
            phiz_path = os.path.join(auto_mesh_working_dir, "phiz.png")
            pylab.savefig(phiz_path)
            pylab.close()
        average_phiz = numpy.mean(array_phiz)
        std_phiz = numpy.std(array_phiz)
        delta_phiz = ny / 2 - average_phiz
    else:
        average_phiz = None
        delta_phiz = None
    # Sample thickness
    list_thickness_index = []
    for omega in [0, 30, 60, 90, 120, 150, 180, 210, 270, 300, 330]:
        str_omega = "%d" % omega
        (array_index, array_upper, array_lower) = dict_loop[str_omega]
        indices = numpy.where((array_index > mesh_xmin) & (array_index < mesh_xmax))
        array_upper = numpy.array(array_upper)[indices]
        array_lower = numpy.array(array_lower)[indices]
        array_thickness = array_upper - array_lower
        # Look for a minima between 100 and 300 pixels from the right:
        array_thickness_crop = array_thickness[-loop_max_width:-loop_min_width]
        indices_thickness_crop = indices[0][-loop_max_width:-loop_min_width]
        if debug:
            pylab.plot(array_upper, "-")
            pylab.plot(array_lower, "*")
            pylab.plot(array_thickness_crop, "+")
            phiz_path = os.path.join(
                auto_mesh_working_dir, "thickness_%s.png" % str_omega
            )
            pylab.savefig(phiz_path)
            pylab.close()
        if len(array_thickness_crop) > 0:
            tmpMin = numpy.argmin(array_thickness_crop)
            #            if tmpMin > 75:
            #                tmpMin = 75
            indexMin = indices_thickness_crop[tmpMin]
            # Minimum mesh length: 150 pixels
            list_thickness_index.append(indexMin)
            logging.debug("Index min: %d for omega = %d" % (indexMin, omega))
    logging.debug("List of thicknesses: %r" % list_thickness_index)
    if not find_largest_mesh and len(list_thickness_index) > 0:
        max_index_thickness = max(list_thickness_index)
        mesh_xmin = max_index_thickness
        logging.debug("Max index for thickess: %d" % max_index_thickness)
    else:
        mesh_xmin = mesh_xmax - loop_max_width
        if mesh_xmin < 0:
            mesh_xmin = 0
    logging.debug("mesh_xmin: %d" % mesh_xmin)
    min_thickness_mesh_ymin = None
    min_thickness_mesh_ymax = None
    max_thickness_mesh_ymin = None
    max_thickness_mesh_ymax = None
    for omega in [0, 30, 60, 90, 120, 150, 180, 210, 270, 300, 330]:
        str_omega = "%d" % omega
        (array_index, array_upper, array_lower) = dict_loop[str_omega]
        indices = numpy.where((array_index > mesh_xmin) & (array_index < mesh_xmax))
        array_upper = numpy.array(array_upper)[indices]
        array_lower = numpy.array(array_lower)[indices]
        if len(array_lower) == 0 or len(array_upper) == 0:
            break
        max_thick = numpy.max(array_upper) - numpy.min(array_lower)
        if max_thick > 400:
            # Ignored
            pass
        else:
            array_index = numpy.arange(mesh_xmin, mesh_xmax - 1)
            if min_thickness is None or min_thickness > max_thick:
                min_thickness = max_thick
                angle_min_thickness = omega
                min_thickness_mesh_ymin = numpy.min(array_lower)
                min_thickness_mesh_ymax = numpy.max(array_upper)
            #             if True:
            #                 pylab.plot(array_upper, '+')
            #                 pylab.plot(array_lower, '+')
            #                 pylab.show()
            if max_thickness is None or max_thickness < max_thick:
                max_thickness = max_thick
                angle_max_thickness = omega
                max_thickness_mesh_ymin = numpy.min(array_lower)
                max_thickness_mesh_ymax = numpy.max(array_upper)
    mesh_xmin -= 50
    logging.debug(
        "Max thickness = %r pxiels at omega %r" % (max_thickness, angle_max_thickness)
    )
    logging.debug(
        "Mesh: xMin=%r, xMax=%r, yMin=%r, yMax = %r"
        % (mesh_xmin, mesh_xmax, max_thickness_mesh_ymin, max_thickness_mesh_ymax)
    )
    logging.debug(
        "Min thickness = %r pixels at omega %r" % (min_thickness, angle_min_thickness)
    )
    logging.debug(
        "Mesh: xMin=%r, xMax=%r, yMin=%r, yMax = %r"
        % (mesh_xmin, mesh_xmax, min_thickness_mesh_ymin, min_thickness_mesh_ymax)
    )
    if mesh_xmin < 10:
        mesh_xmin = 10
    if delta_phiz is not None:
        logging.debug("Average delta phiz: %.4f pixels" % delta_phiz)
        logging.debug("Std delta phiz: %.4f pixels" % std_phiz)
    else:
        logging.warning("Delta phiz could not be determined!")

    x1_pixels = mesh_xmin - nx / 2
    dx_pixels = mesh_xmax - mesh_xmin
    y1_pixels = None
    dy_pixels = None
    angle = None
    if (
        (min_thickness_mesh_ymin is not None)
        and (min_thickness_mesh_ymax is not None)
        and (max_thickness_mesh_ymin is not None)
        and (max_thickness_mesh_ymax is not None)
    ):
        if find_largest_mesh:
            y1_pixels = max_thickness_mesh_ymin - ny / 2
            dy_pixels = max_thickness_mesh_ymax - max_thickness_mesh_ymin
            angle = angle_max_thickness
        else:
            y1_pixels = min_thickness_mesh_ymin - ny / 2
            dy_pixels = min_thickness_mesh_ymax - min_thickness_mesh_ymin
            angle = angle_min_thickness
    debug_message = f"angle={angle}"
    debug_message += f" x1_pixels={x1_pixels}, y1_pixels={y1_pixels}"
    debug_message += f" dx_pixels={dx_pixels}, dy_pixels={dy_pixels}"
    debug_message += f" delta_phiz={delta_phiz}, std_phiz={std_phiz}"
    logging.debug(debug_message)
    return (angle, x1_pixels, y1_pixels, dx_pixels, dy_pixels, delta_phiz, std_phiz)


def plot_img(img, plot_path):
    imgshape = img.shape
    extent = (0, imgshape[1], 0, imgshape[0])
    _ = pyplot.imshow(img, extent=extent)
    pyplot.gray()
    pyplot.colorbar()
    pyplot.savefig(plot_path)
    pyplot.close()
    return


def plotMesh(
    image_path,
    grid_info,
    pixels_per_mm,
    destination_dir,
    sign_phiy=1,
    file_name="snapshot_automesh.png",
    show_plot=False,
):
    (x1_pixels, y1_pixels, dx_pixels, dy_pixels) = gridInfoToPixels(
        grid_info, pixels_per_mm
    )
    img = imageio.imread(image_path, as_gray=True)
    imgshape = img.shape
    extent = (0, imgshape[1], 0, imgshape[0])
    pylab.matshow(img / numpy.max(img), extent=extent)
    pylab.gray()
    ny, nx = imgshape
    if sign_phiy < 0:
        mesh_xmin = nx / 2 - x1_pixels
    else:
        mesh_xmin = nx / 2 + x1_pixels
    mesh_ymin = ny / 2 - y1_pixels
    mesh_xmax = mesh_xmin + dx_pixels
    mesh_ymax = mesh_ymin - dy_pixels
    pylab.plot([mesh_xmin, mesh_xmin], [mesh_ymin, mesh_ymax], color="red", linewidth=2)
    pylab.plot([mesh_xmax, mesh_xmax], [mesh_ymin, mesh_ymax], color="red", linewidth=2)
    pylab.plot([mesh_xmin, mesh_xmax], [mesh_ymin, mesh_ymin], color="red", linewidth=2)
    pylab.plot([mesh_xmin, mesh_xmax], [mesh_ymax, mesh_ymax], color="red", linewidth=2)
    mesh_snap_shot_path = os.path.join(destination_dir, file_name)
    axes = pyplot.gca()
    axes.set_xlim([0, imgshape[1]])
    axes.set_ylim([0, imgshape[0]])
    pylab.savefig(mesh_snap_shot_path, bbox_inches="tight")
    if show_plot:
        pylab.show()
    pylab.close()
    return mesh_snap_shot_path


def gridInfoToPixels(grid_info, pixels_per_mm):
    x1_pixels = grid_info["x1"] * pixels_per_mm
    y1_pixels = grid_info["y1"] * pixels_per_mm
    dx_pixels = grid_info["dx_mm"] * pixels_per_mm
    dy_pixels = grid_info["dy_mm"] * pixels_per_mm
    return (x1_pixels, y1_pixels, dx_pixels, dy_pixels)


def readImage(image_path):
    if image_path.endswith(".png"):
        image = imageio.imread(image_path, as_gray=True)
    elif image_path.endswith(".npy"):
        image = numpy.load(image_path)
    return image


def plotImage(image):
    imgshape = image.shape
    extent = (0, imgshape[1], 0, imgshape[0])
    pyplot.imshow(image, extent=extent)
    pyplot.show()


def plotLoopExam(image, listIndex, listLower, listUpper):
    pyplot.plot(listIndex, listUpper, "+")
    pyplot.plot(listIndex, listLower, "+")
    # imgshape = image.shape
    # extent = (0, imgshape[1], 0, imgshape[0])
    # pyplot.axes(extent)
    pyplot.show()


def cmp(a, b):
    """
    Python3 doesn't have cmp, see:
    https://codegolf.stackexchange.com/questions/49778/how-can-i-use-cmpa-b-with-python3
    """
    return (a > b) - (a < b)
