import matplotlib.pyplot as plt
import numpy as np
import cv2 as cv
import os
import sys
from astropy.io import fits


class FindLensTemplate:

    def __init__(self, image_dir, image_prefix, temp_path,
                 matchTemplate=cv.matchTemplate,
                 method=cv.TM_CCOEFF_NORMED,
                 isMax=True, rescale="log",
                 ndtype="float32",
                 cutout = False,
                 y_slice = None,
                 x_slice = None):
        """
        This module is to find similar patterns between different band of same lensing fields
        :param image_dir: str, the directory of image lists
        :param temp_name: str, the name of template or template image in the same directory of images
        :param matchTemplate: callable, matching function
        :param method: None, the parameters of function
        :param isMax: Bool, output the maximum or minimum result of matching, default is True
        :param rescale: str, mathematical functions from numpy to rescale the images.
        :param ndtype: str, the dtype of image and template
        """
        self.image_dir = image_dir
        self.image_prefix = image_prefix
        self.temp_path = temp_path
        self.method = method
        self.isMax = isMax
        self.rescale = rescale
        self.matchTemplate = matchTemplate
        self.ndtype = ndtype
        self.bands_order = ['z', 'i', 'r', 'g']

        self.cutout = cutout #[1098:1218,544:664]
        self.y_slice = y_slice
        self.x_slice = x_slice

        if rescale == "linear":
            self.rescaler = eval("np.array")
        elif rescale == "power":
            self.rescaler = eval("np.square")
        else:
            try:
                self.rescaler = eval("np." + self.rescale)
            except:
                print("scale is not a numpy function.")



    def MethodMatch(self, image_path, temp):

        """
        Find the most similar temp pattern in single image
        :param image_path: str, the path of target image
        :param temp: str or ndarray, the name of template image, if ndarray, it is the image of template
        :return: ndarray, the most likely position in the image
        """

        # open the image fits file and rescale/normalize the image
        image_file = fits.open(image_path)   # open fits file
        image = image_file[0].data      # load fits data
        if image is None:
            image = image_file[1].data #in case the image data is stored in [1]
        image = image.astype(self.ndtype)
        image_min = image.min()
        #preparing the data for algorithms that require non-negative values or are sensitive to values being too close to zero.
        image = image - image_min * (1 - np.sign(image_min) * 0.001)    # parallel the image values
        image = self.rescaler(image)    # rescale the image
        image_min = image.min()
        image_max = image.max()
        image = (image - image_min) / (image_max - image_min) * 255 # #to ensure to normalize the images to a range of [0,255]
        # image = (image) / (image_max) * 255  #
        img = image.copy()

        if self.cutout:
            img = img[self.y_slice, self.x_slice]

        # checking the type of template and matching the image
        if isinstance(temp, str):
            template = fits.open(temp)[0].data
            res = self.matchTemplate(img, template, method=self.method)

        elif isinstance(temp, np.ndarray):

            res = self.matchTemplate(img, temp, method=self.method)

        # extract the minimum and maximum piont of residual image
        min_val, max_val, min_loc, max_loc = cv.minMaxLoc(res)

        if self.isMax:
            if self.cutout:
                return (max_loc[0] + self.x_slice.start, max_loc[1]+self.y_slice.start)
            else:
                return max_loc
        else:
            return min_loc

    def SelectFiles(self, file_list, prefix=2, seperator="_"):
        #add, fix the sequence z,i,r,g
        selected_file_list = []
        for file in file_list:
            file_split = file.split(".")[0].split(seperator)
            try:
                band = file_split.pop(prefix)
                select_file_name = seperator.join(file_split)

                if select_file_name == self.image_prefix:
                    selected_file_list.append((file, band))
            except:
                pass
        selected_file_list.sort(key=lambda x: self.bands_order.index(x[1])) #x[1] extracts the band part of the tuple; using sort() to Elements are ordered so that their bands appear in the same order as they do in self.bands_order
        # Return only the file names, excluding the band identifiers
        return [file for file, band in selected_file_list]

    def MulFileMatch(self, progress=10, prefix=2, seperator="_"):
        """

        :param progress:
        :param prefix:
        :param seperator:
        :return:
        """

        template = fits.open(self.temp_path)[0].data
        template = template.astype(self.ndtype)

        position_dict = {}

        file_list = os.listdir(self.image_dir)
        file_list = self.SelectFiles(file_list, prefix=prefix, seperator=seperator) #z,i,r,g sequence
        dir_len = len(file_list)

        if progress!=0:
            for file_n in range(dir_len):

                file_name = file_list[file_n]

                file_path = os.path.join(self.image_dir, file_name)

                position = self.MethodMatch(file_path, template)

                position_dict[file_list[file_n]] = position

                sys.stdout.write('\r')
                # the exact output you're looking for:
                sys.stdout.write(("[%-"+str(dir_len*progress)+"s] %d%% %d/%d")\
                                 % ('='*(file_n+1)*progress, (file_n+1)/(dir_len)*100, file_n+1, dir_len))
                sys.stdout.flush()

        elif progress == 0:
            for file_n in range(dir_len):
                file_name = file_list[file_n]
                if (file_name != self.temp_name):
                    file_path = os.path.join(self.image_dir, file_name)

                    position = self.MethodMatch(file_path, template)

                    position_dict[file_list[file_n]] = position
        else:
            print("Progress only accepts an integer.")

        return position_dict

    def MulMethodMatch(self, image, temp, methods=None, isMax=True):
        """

        :param image:
        :param temp:
        :param methods:
        :param isMax:
        :return:
        """

        if methods is None:
            methods = [cv.TM_CCOEFF_NORMED]

        res_ms = []

        img = image.copy()

        for method in methods:
            res_m = cv.matchTemplate(img, temp, method=method)

            res_ms.append(res_m)

        res_ms = np.array(res_ms)

        min_val, max_val, min_loc, max_loc = cv.minMaxLoc(res_ms)

        if isMax:
            return max_loc
        else:
            return min_loc
