from astropy.io import fits

class SubregionSaver:
    def __init__(self, prefix):
        self.prefix = prefix

    def save_subregion_to_fits(self, bands_data, output_path,
                               suffix = 'cut', #400, 800
                               save_bool=False):
        #bands_data
        for band, data in bands_data.items():
            hdu = fits.PrimaryHDU(data)
            hdul = fits.HDUList([hdu])

            if save_bool:
                # Save the file with '_cut.fits' suffix
                filename_cut = f'{output_path}/{self.prefix}_{band}_{suffix}.fits'
                hdul.writeto(filename_cut, overwrite=True)



# # Example usage
# output_path = '../../lenstronomy/dataset/data_prepared_to_reconstruct'
# prefix = 'coadd_DESJ0010-4315'
# bands_data = {
#     'z': Image_lens[0],  # Assuming Image_lens is previously defined
#     'i': Image_lens[1],
#     'r': Image_lens[2],
#     'g': Image_lens[3],
# }
#
# # Creating an instance of the class with a specific prefix
# saver = SubregionSaver(prefix)
#
# # Using the method with the new parameters
# saver.save_subregion_to_fits(bands_data, output_path, save_as_cut=True)
