"""
Utility methods for cutting image stamps from fits files. Various methods for
specifying which region to cut.
"""
from __future__ import division
import astropy.io.fits as fits
import astropy.wcs as wcs
import pyregion
import warnings
from numpy import where, isfinite, array
import re

_stamp_slice_key = 'STAMPSLC'
_stamp_slice_comment = 'Numpy array slice [Y,X] of stamp area in original image'


def _slice_to_string(nd_slice):
    """
    Convert an n-dimensional numpy slice into a string. Use only start and stop
    (step not used).
    :param nd_slice: Sequence of slice objects (e.g. for numpy indexing)
    :rtype: String representation of slice: [a:b, ... , y:z]
    """
    slice_str = ['{0:d}:{1:d}'.format(nd_slice[dim].start, nd_slice[dim].stop)
                 for dim in xrange(len(nd_slice))]
    return '[' + ','.join(slice_str) + ']'


def _string_to_slice(slice_str):
    """
    Convert string specifying a numpy n-dimensional slice into
    :param slice_str: String in usual numpy format: [a:b, ... , y:z]
    :rtype: Tuple of slices represented in slice_string
    """
    slice_str = slice_str.strip('[]')
    slice_str = re.split('\D+', slice_str)
    slice_str = [int(num) for num in slice_str]
    return tuple([slice(slice_str[2*dim], slice_str[2*dim + 1])
                  for dim in xrange(len(slice_str) // 2)])


def _stamp_from_regions(hdu_list, regions, extension=0, **kwargs):
    """
    Cut a stamp using a list of pyregion regions to generate the mask

    :param hdu_list: Source image (FITS HDUList)
    :param regions: List of pyregion Shapes, or pyregion ShapeList
    :param extension: extension of the HDU list to cut data from
    :param kwargs: Additional arguments passed to _stamp_from_mask
    :rtype: FITS PrimaryHDU with the original image header information.
    """
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
    extents = tuple([slice(masked_indexes[dim].min(),
                           masked_indexes[dim].max()+1)
                     for dim in xrange(len(masked_indexes))])
    mask = mask[extents]

    # copy the data of the hduList
    new_data = hdu_list[extension].data[extents].copy()
    # Replace masked-out areas with masked_value
    if mask_inf_nan:
        mask &= isfinite(new_data)
    if masked_value is not None:
        new_data[~mask] = masked_value
    new_hdu = fits.PrimaryHDU(header=hdu_list[extension].header, data=new_data)

    wcs_file = wcs.WCS(hdu_list[extension].header)
    coord_types = [axis['coordinate_type'] for axis
                   in wcs_file.get_axis_types()]
    if all([ctype == 'celestial' for ctype in coord_types]):
        # Re-Center WCS if celestial coordinate transform is present
        # Center is average of min and max pixels. Fits standard is 1-based and
        # X,Y. Extents slice is numpy 0-based, stop non-inclusive, and Y,X
        center = array((0.5*(extents[1].start + 1 + extents[1].stop),
                        0.5*(extents[0].start + 1 + extents[0].stop)))
        # The second argument here specifies the coordinate origin (1 for FITS)
        wcs_file.wcs.crval = wcs_file.all_pix2world([center], 1)[0]
        # Center in new image. Equivalently, 0.5*(reversed(data.shape)+1)
        wcs_file.wcs.crpix = center - (extents[1].start, extents[0].start)
        new_hdu.header.extend(wcs_file.to_header(), update=True)
    elif all([ctype is None for ctype in coord_types]):
        # When no WCS information is present, skip coordinate re-centering
        # without comment to the user.
        pass
    else:
        warnings.warn('File {} contains a mixture of celestial and non-'
                      'celestial WCS axes. No WCS re-centering will be '
                      'performed.'.format(hdu_list.filename()),
                      wcs.AstropyUserWarning)
    new_hdu.update_header()

    orig_sec_str = _slice_to_string(extents)
    new_hdu.header[_stamp_slice_key] = (orig_sec_str, _stamp_slice_comment)
    return new_hdu


def _cut_stamps_regions(source_file, region_file, source_ext=0, file_prefix='',
                        selected_labels=None, **kwargs):
    """
    Cuts stamps from a FITS image based on a DS9 region file. See docstring for
    cut_stamps for argument explanations
    """
    # Open region file first, so parse errors will fail before trying to open
    # a potentially large FITS file
    regions = pyregion.open(region_file)

    # Open cutting image
    hdu_list = fits.open(source_file, mode='readonly', ignore_missing_end=True)

    # Translate to image coordinates, since we're working in pixel space
    regions = regions.as_imagecoord(hdu_list[source_ext].header)

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
                                    extension=source_ext, **kwargs)
        filename = file_prefix + name + '.fits'
        stamp.writeto(filename, clobber=True)
        all_stamp_names += [filename]
    hdu_list.close()

    return all_stamp_names


def _cut_stamps_labels(source_file, label_file, source_ext=0, file_prefix='',
                       selected_labels=None, **kwargs):
    """
    Cuts stamps from a FITS image based on a FITS-format label image (e.g.
    sextractor segmentation map). See docstring for cut_stamps for argument
    explanations
    """
    label_data = fits.getdata(label_file)
    # If no specific labels were requested, cut stamps for every object
    if selected_labels is None:
        selected_labels = xrange(1, label_data.max()+1)

    # Open cutting image
    hdu_list = fits.open(source_file, mode='readonly', ignore_missing_end=True)

    all_stamp_names = []
    for obj_num in selected_labels:
        name = str(obj_num)
        stamp = _stamp_from_label(hdu_list, label_data, label_num=obj_num,
                                  extension=source_ext, **kwargs)
        filename = file_prefix + name + '.fits'
        stamp.writeto(filename, clobber=True)
        all_stamp_names += [filename]
    hdu_list.close()

    return all_stamp_names


def cut_stamps(source_file, region_file, region_format='region', source_ext=0,
               file_prefix='', selected_labels=None, masked_value=None,
               mask_inf_nan=False):
    """
    Cuts stamps from imgfile for each region specified in a regfile. Two region
    formats are currently supported:
    "region": Specifies regfile is in SAOImage DS9 region format. Regions with
        the same name will be combined into a single mask. This can be used to
        create composite, multi-component regions.
    "label": Specifies regfile is a label image where each object's pixels are
        flagged by a unique non-zero integer (e.g., sextractor segmentation map)
        Datacubes with more than 2 dimensions are supported for label images.

    :param source_file: Filename of the image to cut from
    :param region_file: Filename specifying the regions to cut. Expected format
        is controlled by the region_format argument.
    :param region_format: The file format of regfile. Current options:
        "region", "ds9" - for SAOImage DS9 region files
        "label", "segmentation", "sextractor" - for FITS images with integer
        unique labels for each object (i.e. sextractor segmentation maps)
    :param source_ext: FITS extension that contains the image data to cut
    :param file_prefix: Base filename, to which region name or number is
        appended. E.g., if cutting from weight maps, perhaps "wht_"
    :param selected_labels: List of specific label numbers (or DS9 region text
        attributes) to output. If None, stamps are cut for all labels/regions.
    :param masked_value: Value with which to fill masked-out regions of stamps.
        E.g. when cutting weight map stamps, masked_value=0 may be useful.
    :param mask_inf_nan: Whether to replace inf/nan pixels with masked_value
    :rtype: List of filenames of created stamps.
    """
    region_format = region_format.lower()
    if region_format in ('region', 'ds9'):
        cut_function = _cut_stamps_regions
    elif region_format in ('label', 'segmentation', 'sextractor'):
        cut_function = _cut_stamps_labels
    else:
        raise ValueError('Unknown region file format')

    return cut_function(source_file, region_file, source_ext=source_ext,
                        file_prefix=file_prefix,
                        selected_labels=selected_labels,
                        masked_value=masked_value,
                        mask_inf_nan=mask_inf_nan)


def paste_stamps(target_file, stamp_files, target_ext=0, stamps_ext=0,
                 output_name=None):
    """
    Pastes a collection of cut stamp images back into an original image. For
    instance, if you cut a stamp, model out a galaxy or star using e.g. galfit,
    and want to propagate the modified stamp back into the original image.

    :param target_file: Image that the stamps will be pasted into
    :param stamp_files: Stamp images that will be pasted into dest_file. Can
        either be a single filename or a list of filenames. Must retain the
        STAMPSLC FITS keyword originally added during cutting.
    :param target_ext: FITS extension of dest_file to paste into
    :param stamps_ext: FITS extension of the stamp images to copy from
    :param output_name: Name of the new file with pasted stamps. If None,
        the suffix "_pasted" will be appended to dest_file. E.g. original.fits
        becomes original_pasted.fits. If output_name is the same as dest_file,
        dest_file will be overwritten.
    :rtype: Filename of the modified copy of target_file
    """
    if isinstance(stamp_files, basestring):
        stamp_files = [stamp_files]
    if output_name is None:
        output_name = target_file.replace('.fits', '_pasted.fits')
    destination = fits.open(target_file)
    for stamp_file in stamp_files:
        stamp = fits.open(stamp_file)
        paste_slice = stamp[stamps_ext].header[_stamp_slice_key]
        paste_slice = _string_to_slice(paste_slice)
        slice_size = tuple([sl.stop - sl.start for sl in paste_slice])
        if slice_size != stamp[stamps_ext].data.shape:
            raise ValueError('Shape of stamp image does not match shape of '
                             '{0} header keyword.'.format(_stamp_slice_key))
        destination[target_ext].data[paste_slice] = stamp[stamps_ext].data
        stamp.close()
    destination.writeto(output_name, clobber=True)
    destination.close()

    return output_name
