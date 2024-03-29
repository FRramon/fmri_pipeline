# first line: 687
def resample_to_img(
    source_img,
    target_img,
    interpolation="continuous",
    copy=True,
    order="F",
    clip=False,
    fill_value=0,
    force_resample=False,
):
    """Resample a Niimg-like source image on a target Niimg-like image.

    No registration is performed: the image should already be aligned.

    .. versionadded:: 0.2.4

    Parameters
    ----------
    source_img : Niimg-like object
        See :ref:`extracting_data`.
        Image(s) to resample.

    target_img : Niimg-like object
        See :ref:`extracting_data`.
        Reference image taken for resampling.

    interpolation : str, default='continuous'
        Can be 'continuous', 'linear', or 'nearest'. Indicates the resample
        method.

    copy : bool, default=True
        If True, guarantees that output array has no memory in common with
        input array.
        In all cases, input images are never modified by this function.

    order : "F" or "C", default="F"
        Data ordering in output array. This function is slightly faster with
        Fortran ordering.

    clip : bool, default=False
        If False (default) no clip is performed.
        If True all resampled image values above max(img)
        and under min(img) are cllipped to min(img) and max(img).

    fill_value : float, default=0
        Use a fill value for points outside of input volume.

    force_resample : bool, default=False
        Intended for testing, this prevents the use of a padding optimization.

    Returns
    -------
    resampled: nibabel.Nifti1Image
        input image, resampled to have respectively target image shape and
        affine as shape and affine.

    See Also
    --------
    nilearn.image.resample_img

    """
    target = _utils.check_niimg(target_img)
    target_shape = target.shape

    # When target shape is greater than 3, we reduce to 3, to be compatible
    # with underlying call to resample_img
    if len(target_shape) > 3:
        target_shape = target.shape[:3]

    return resample_img(
        source_img,
        target_affine=target.affine,
        target_shape=target_shape,
        interpolation=interpolation,
        copy=copy,
        order=order,
        clip=clip,
        fill_value=fill_value,
        force_resample=force_resample,
    )
