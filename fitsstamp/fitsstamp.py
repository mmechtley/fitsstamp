"""
Utility methods for cutting image stamps from fits files. Various methods for
specifying which region to cut.
"""
import os
import astropy.io.fits as fits
import astropy.wcs as wcs
import pyregion
from numpy import where, isfinite, array

_stampsec_key = 'STAMPSEC'
_stampsec_comment = 'Section of original image from which stamp was extracted'


def _stamp_from_regions(hdu_list, regions, extension=0, masked_value=None,
                        mask_inf_nan=False):
    """
    Cut a stamp using mask data (such as from pyregion spatial filter)
    :param hdu_list: source image (pyfits HDU list)
    :param regions: list of pyregion Shapes, or pyregion ShapeList
    :param extension: extension of the HDU list to cut data from
    :param masked_value: The value to replace masked-out values with.
    :param mask_inf_nan: Replace inf and nan values with masked_value
    :return: Pyfits HDU list with the original header information.
            You must save this HDUList to a file yourself!
    """
    # TODO: Something in this function (or writeto?) is slower than expected
    data_shape = hdu_list[extension].data.shape
    regions = pyregion.ShapeList(regions)
    mask = regions.get_filter().mask(data_shape)

    masked_indexes = where(mask)
    extents = (slice(masked_indexes[0].min(), masked_indexes[0].max()+1),
               slice(masked_indexes[1].min(), masked_indexes[1].max()+1))
    mask = mask[extents]

    # copy the data of the hduList
    new_data = hdu_list[extension].data[extents].copy()
    # Replace masked-out areas with masked_value
    if mask_inf_nan:
        mask &= isfinite(new_data)
    if masked_value is not None:
        new_data[~mask] = masked_value
    new_hdu = fits.PrimaryHDU(header=hdu_list[extension].header, data=new_data)

    # Re-Center WCS
    wcs_file = wcs.WCS(hdu_list[extension].header)
    # center is average of min and max pixels. Fits standard is 1-based
    center = array((0.5*(extents[1].start + extents[1].stop - 1) + 1,
                    0.5*(extents[0].start + extents[0].stop - 1) + 1))
    wcs_file.wcs.crval = wcs_file.all_pix2sky([center], 1)[0]
    # extents are 0-based, so no need to subtract 1 from bottom/left pixels
    wcs_file.wcs.crpix = center - (extents[1].start, extents[0].start)
    new_hdu.header.extend(wcs_file.to_header(), update=True)
    new_hdu.update_header()

    orig_sec_str = '[{0:d}:{1:d},{2:d}:{3:d}]'.format(
        extents[1].start, extents[1].stop, extents[0].start, extents[0].stop)
    new_hdu.header.update(_stampsec_key, orig_sec_str,
                          comment=_stampsec_comment)
    return new_hdu


# TODO: Cut stamps using label image (e.g. source extractor segmentation map)
def cut_regions_from_image(imgfile, regfile, extension=0, base_name='',
                           masked_value=None, mask_inf_nan=False):
    """
    Cuts stamps for each region specified in a ds9-format region file. Regions
    with the same name will be combined into a single mask. Use this to create
    composite, multi-component regions.
    :param imgfile: Filename of the image to cut from
    :param regfile: ds9 region file to get the cut points
    :param extension: extension that contains the image data to cut
    :param base_name: base filename to which region name or number is appended
    :param masked_value: value with which to fill masked-out regions of stamps
    :param mask_inf_nan: Whether to convert inf/nan pixels to masked_value
    """
    # Only accept local files for cutting -- none of this weird URL stuff
    if not os.path.exists(imgfile) or not os.path.exists(regfile):
        raise IOError('Image or Region file does not exist. Perhaps you '
                      'mistyped it?')
    
    # Open region file first, to allow malformed regions to fail before
    # opening a giant image
    regions = pyregion.open(regfile)

    # Open cutting image
    hdu_list = fits.open(imgfile, mode='readonly', ignore_missing_end=True)

    # Translate to image coordinates
    regions = regions.as_imagecoord(hdu_list[extension].header)

    # If the region had a text label, use that as the default file name
    names = [reg.attr[1].get('text', str(i)) for i, reg in enumerate(regions)]
    
    all_stamps = []
    for name in set(names):
        name_regs = [reg for reg, nm in zip(regions, names) if nm == name]
        stamp = _stamp_from_regions(hdu_list, name_regs, extension=extension,
                                    masked_value=masked_value,
                                    mask_inf_nan=mask_inf_nan)
        filename = base_name + name + '.fits'
        stamp.writeto(filename, clobber=True)
        all_stamps += [filename]
    
    hdu_list.close()

    return all_stamps
