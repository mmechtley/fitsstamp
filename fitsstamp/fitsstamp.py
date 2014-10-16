"""
Utility methods for cutting image stamps from fits files. Various methods for
specifying which region to cut.
"""
import astropy.io.fits as fits
import astropy.wcs as wcs
from numpy import where, isfinite, array
try:
    import pyregion
except ImportError:
    pyregion = None

_stampsec_key = 'STAMPSEC'
_stampsec_comment = 'Section of original image from which stamp was extracted'


def _stamp_from_regions(hdu_list, regions, extension=0, **kwargs):
    """
    Cut a stamp using a list of pyregion regions to generate the mask

    :param hdu_list: Source image (FITS HDUList)
    :param regions: List of pyregion Shapes, or pyregion ShapeList
    :param extension: extension of the HDU list to cut data from
    :param kwargs: Additional arguments passed to _stamp_from_mask
    :rtype: FITS PrimaryHDU with the original image header information.
    """
    # TODO: Something in this function (or writeto?) is slower than expected
    data_shape = hdu_list[extension].data.shape
    regions = pyregion.ShapeList(regions)
    mask = regions.get_filter().mask(data_shape)
    return _stamp_from_mask(hdu_list, mask, extension=extension, **kwargs)


def _stamp_from_label(hdu_list, label_image, label_num=1, **kwargs):
    """
    Cut stamp using a label image (like sextractor segmentation map), finding
    those pixels that belong to a certain label number.

    :param hdu_list: Source image (FITS HDUList)
    :param label_image: Numpy array with same shape as hdu_list image, where
        each object's pixels are marked by a unique integer
    :param label_num: Integer label of the object to be selected
    :param kwargs: Additional arguments passed to _stamp_from_mask
    :rtype: FITS PrimaryHDU with the original image header information.
    """
    mask = label_image == label_num
    return _stamp_from_mask(hdu_list, mask, **kwargs)


def _stamp_from_mask(hdu_list, mask, extension=0, masked_value=None,
                     mask_inf_nan=False):
    """
    Cut a stamp using a numpy mask array (True for selected pixels, False for
    ignored pixels)

    :param hdu_list: Source image (FITS HDUList)
    :param mask: Numpy mask array (bool, same shape as image, selected pixels
        are True and non-selected pixels are False)
    :param extension: Extension of the HDUList to cut data from
    :param masked_value: The value to replace non-selected pixels with.
    :param mask_inf_nan: Replace inf and nan values with masked_value
    :rtype: FITS PrimaryHDU with the original image header information.
    """
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
    # center is average of min and max pixels. Fits standard is 1-based and X,Y
    # whereas extents slice is numpy 0-based and Y,X
    center = array((0.5*(extents[1].start + 1 + extents[1].stop),
                    0.5*(extents[0].start + 1 + extents[0].stop)))
    # The second argument here specifies the coordinate origin (1 for FITS)
    wcs_file.wcs.crval = wcs_file.all_pix2world([center], 1)[0]
    # Center in new image. Alternatively, 0.5*(reversed(shape)+1)
    wcs_file.wcs.crpix = center - (extents[1].start, extents[0].start)
    new_hdu.header.extend(wcs_file.to_header(), update=True)
    new_hdu.update_header()

    orig_sec_str = '[{0:d}:{1:d},{2:d}:{3:d}]'.format(
        extents[1].start, extents[1].stop, extents[0].start, extents[0].stop)
    new_hdu.header.update(_stampsec_key, orig_sec_str,
                          comment=_stampsec_comment)
    return new_hdu


def _cut_stamps_regions(imgfile, regfile, extension=0, base_name='',
                        selected_labels=None, **kwargs):
    """
    Cuts stamps from a FITS image based on a DS9 region file. See docstring for
    cut_stamps for argument explanations
    """
    # Open region file first, so parse errors will fail before trying to open
    # a potentially large FITS file
    regions = pyregion.open(regfile)

    # Open cutting image
    hdu_list = fits.open(imgfile, mode='readonly', ignore_missing_end=True)

    # Translate to image coordinates, since we're working in pixel space
    regions = regions.as_imagecoord(hdu_list[extension].header)

    # If a region has a text label, use that as the default file name
    # reg.attr[1] is a dictionary of extended region attributes
    names = [reg.attr[1].get('text', str(num))
             for num, reg in enumerate(regions)]
    if selected_labels is None:
        selected_labels = set(names)

    all_stamp_names = []
    for name in set(names) & selected_labels:
        regs_for_name = [reg for reg, nm in zip(regions, names) if nm == name]
        stamp = _stamp_from_regions(hdu_list, regs_for_name,
                                    extension=extension, **kwargs)
        filename = base_name + name + '.fits'
        stamp.writeto(filename, clobber=True)
        all_stamp_names += [filename]
    hdu_list.close()

    return all_stamp_names


def _cut_stamps_labels(imgfile, labelfile, extension=0, base_name='',
                       selected_labels=None, **kwargs):
    """
    Cuts stamps from a FITS image based on a FITS-format label image (e.g.
    sextractor segmentation map). See docstring for cut_stamps for argument
    explanations
    """
    label_data = fits.getdata(labelfile)
    # If no specific labels were requested, cut stamps for every object
    if selected_labels is None:
        selected_labels = xrange(1, label_data.max()+1)

    # Open cutting image
    hdu_list = fits.open(imgfile, mode='readonly', ignore_missing_end=True)

    all_stamp_names = []
    for obj_num in selected_labels:
        name = str(obj_num)
        stamp = _stamp_from_label(hdu_list, label_data, label_num=obj_num,
                                  extension=extension, **kwargs)
        filename = base_name + name + '.fits'
        stamp.writeto(filename, clobber=True)
        all_stamp_names += [filename]
    hdu_list.close()

    return all_stamp_names


def cut_stamps(imgfile, regfile, region_format='region', extension=0,
               base_filename='', selected_labels=None, masked_value=None,
               mask_inf_nan=False):
    """
    Cuts stamps from imgfile for each region specified in a regfile. Two region
    formats are currently supported:
    "region": Specifies regfile is in SAOImage DS9 region format. Regions with
        the same name will be combined into a single mask. This can be used to
        create composite, multi-component regions.
    "label": Specifies regfile is a label image where each object's pixels are
        flagged by a unique non-zero integer (e.g., sextractor segmentation map)

    :param imgfile: Filename of the image to cut from
    :param regfile: Filename specifying the regions to cut. Expected format is
        controlled by the region_format argument.
    :param region_format: The file format of regfile. Current options:
        "region" - for SAOImage DS9 region files
        "label" - for FITS images with integer unique labels for each object
            (i.e. sextractor segmentation maps)
    :param extension: FITS extension that contains the image data to cut
    :param base_filename: Base filename, region name or number is appended.
        E.g., if cutting science or weight maps, perhaps "sci_" or "wht_"
    :param selected_labels: List of specific label numbers (or DS9 region text
        attributes) to output. Otherwise stamps are cut for all labels/regions.
    :param masked_value: Value with which to fill masked-out regions of stamps.
        E.g. when cutting weight map stamps, masked_value=0 may be useful.
    :param mask_inf_nan: Whether to replace inf/nan pixels with masked_value
    :rtype: List of filenames of created stamps.
    """
    if region_format in ('region', 'DS9'):
        cut_function = _cut_stamps_regions
    elif region_format in ('label', 'segmentation', 'sextractor'):
        cut_function = _cut_stamps_labels
    else:
        raise ValueError('Unknown region file mode')

    return cut_function(imgfile, regfile, extension=extension,
                        base_name=base_filename,
                        selected_labels=selected_labels,
                        masked_value=masked_value,
                        mask_inf_nan=mask_inf_nan)