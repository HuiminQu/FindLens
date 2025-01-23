import numpy as np
import cv2 as cv
import os
from astropy.io import fits
from matplotlib import pyplot as plt
from FindLensTemplate import FindLensTemplate
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

class VisualTemp:
    def __init__(self, data_dir, position_dict, ncols=2, nrows=2, w=50, h=50, **kwargs):

        self.data_dir = data_dir
        self.position_dict = position_dict
        self.len = len(position_dict)
        self.w = w
        self.h = h
        self.ncols = ncols
        self.nrows = nrows

        if self.len > (self.ncols*self.nrows):
            print("Error, The lens of position dictionary is out of range of subplots.")
            exit()

    def ShowSubplots(self, cmap="rainbow", figsize=(8, 8)):
        #wavelength from long to short: z,i,r,g
        file_name_list = list(self.position_dict)
        position_list = list(self.position_dict.values())

        fig, ax = plt.subplots(self.nrows, self.ncols, figsize=figsize)

        for r in range(self.nrows):
            for c in range(self.ncols):
                try:
                    image = fits.open(os.path.join(self.data_dir,
                                file_name_list[c+self.ncols*r]))
                    image = image[0].data
                    if image is None:
                        image = image[1].data  # in case the image data is stored in [1]
                    position = position_list[c+self.ncols*r]
                    img = ax[r, c].imshow(image[position[1]:position[1]+self.h, #y
                                    position[0]:position[0]+self.w], origin='lower', cmap=cmap)
                    ax[r, c].set_title(file_name_list[c+self.ncols*r])
                    fig.colorbar(img)
                except:
                    ax[r,c].axis("off")
        fig.show()

        return fig, ax


    def ShowPosition(self, cmap="rainbow", dotcolor="red", maker="*", size=30, figsize=(6, 8)):

        file_name_list = list(self.position_dict)
        position_list = list(self.position_dict.values())

        fig, ax = plt.subplots(self.nrows, self.ncols, figsize=figsize)

        for r in range(self.nrows):
            for c in range(self.ncols):
                try:
                    image = fits.open(os.path.join(self.data_dir,
                                file_name_list[c+self.ncols*r]))
                    image = image[0].data
                    if image is None:
                        image = image[1].data  # in case the image data is stored in [1]
                    position = position_list[c+self.ncols*r]
                    img = ax[r, c].imshow(image, origin='lower', cmap=cmap)
                    ax[r, c].scatter(position[0], position[1],
                                     c=dotcolor, marker=maker, s=size)
                    ax[r, c].set_title(file_name_list[c+self.ncols*r])

                except:
                    ax[r, c].axis("off")
        fig.show()

        return fig, ax

    def lensSubtractbkg(self):
        Image_lens = [[], [], [], []]
        Data_bkgsub = [[], [], [], []]

        file_name_list = list(self.position_dict)
        position_list = list(self.position_dict.values())

        for r in range(self.nrows):
            for c in range(self.ncols):
                try:
                    image = fits.open(os.path.join(self.data_dir,
                                file_name_list[c+self.ncols*r]))
                    image = image[0].data
                    if image is None:
                        image = image[1].data  # in case the image data is stored in [1]
                    
                    sigma_clip = SigmaClip(sigma=3.)
                    coverage_mask = (image == 0)
                    bkg_estimator = MedianBackground()
                    bkg = Background2D(image, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip,
                                       bkg_estimator=bkg_estimator, coverage_mask=coverage_mask, fill_value=0.0)
                    data_bkgsub = image.copy()
                    data_bkgsub = data_bkgsub - bkg.background
                    Data_bkgsub[c + self.ncols * r] = data_bkgsub
                    
                    position = position_list[c+self.ncols*r]
                    Image_lens[c + self.ncols * r] = data_bkgsub[position[1]:position[1] + self.h,
                                                     position[0]:position[0] + self.w] #y,x; while position gives in sequence [x,y]

                except:
                    pass
        self.Data_bkgsub = Data_bkgsub
        return Image_lens

    def Subtractbkg_arbitrary_size(self, size=400):
        Image_arbitrary_size = [[], [], [], []]
        position_list = list(self.position_dict.values())
        for r in range(self.nrows):
            for c in range(self.ncols):
                position = position_list[c + self.ncols * r]
                data_bkgsub = self.Data_bkgsub[c + self.ncols * r]
                Image_arbitrary_size[c + self.ncols * r] = data_bkgsub[int(position[1]-(size-self.h)/2):int(position[1]+self.h+(size-self.h)/2),
                                                           int(position[0]-(size-self.w)/2):int(position[0]+self.w+(size-self.w)/2)]
        return Image_arbitrary_size