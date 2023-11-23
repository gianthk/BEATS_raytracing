# -*- coding: utf-8 -*-
"""
This module provides a few wrapper functions around lmfit to easily fit
profiles to a selection of given functions and optional background functions.

Profiles can also be directly extracted from image data, either along the
horizontal or vertical image directions.

Convenience functions for directly grabing a new image and running the fitting
on this image are provided.

Note
----
The wrapping of lmfit into these easy to use functions does shadow a lot of its
advanced functionality. Initial parameter guesses are determined automatically
and may not work in all cases. If fitting with the functions provided in this
module fails, it is recommended to fall back to the (also relatively easy to
use) native functionality of lmfit.

"""


import matplotlib.pyplot as plt
import numpy as np

from lmfit.models import (ConstantModel, LinearModel, GaussianModel,
    LorentzianModel, PseudoVoigtModel, RectangleModel, StepModel)

# from tomoalign.scan import grab_image
# from tomoalign.utils import get_beamline_config, get_camera


def get_profile(im, axis=0, width=None, position=None, trim=0):
    """
    Extract a horizontal or vertical linear profile from image data

    Parameters
    ----------
    im : 2D numpy array-like
        The image data which is to be used in the fit.
    axis : str or int, optional
        The axis direction along which the profile is to be computed. Valid
        options are:

        * 'v', 'vertical', or 0: vertical profile
        * 'h', 'horizontal', or 1: horizontal profile

        (default = 0 --> vertical)
    width : int, optional
        Full width of the integration window perpendicular to the profile
        direction in number of pixels. If `None` is specified, the width is
        taken to be the full image size. The mean value over the window width
        is calculated. (default = None)
    position : int, optional
        The center position of the profile in pixel coordinates. If `None` is
        specified, the position is set to the middle of the image.
        (default = None)
    trim : int or [int, int], optional
        The number of pixels to trim from the image edge at the start and end
        of the profile. If only one number is specified, the profile is
        trimmed symmetrically. (default = 0)

    Returns
    -------
    x : 1D array-like
        The independent variable (pixel position) of the profile data.
    y : 1D array-like
        The dependent variable (intensity value) data of the profile data.

    """

    if type(axis) is str:
        if axis.startswith('h'):
            axis = 1
        else:
            axis = 0

    if axis > 0:
        im = im.T

    if position is None:
        position = int(im.shape[1] / 2.0)
    if width is None:
        width = im.shape[1]
    w = int(width / 2.0)

    if not type(trim) is list:
        trim = list([trim])
    if not len(trim) == 2:
        trim = [trim[0], trim[0]]
    if np.any(trim):
        if trim[1] > 0:
            im = im[trim[0]:-trim[1], :]
        else:
            im = im[trim[0]:, :]

    profile = np.mean(im[:, position-w+1:position+w], axis=1)
    x = np.arange(len(profile)) + trim[0]

    return x, profile


def fit_beam_edge(direction='vertical', image=None, width=40, trim=0):
    """
    Fit a step-edge beam profile on a constant background to the central 40
    pixel wide region of the image along the specified direction.

    Plot the result and display the fitted center position.

    Parameters
    ----------
    direction : str or int
        The axis direction along which the profile is to be computed. Valid
        options are:

        * 'v', 'vertical', or 0: vertical profile
        * 'h', 'horizontal', or 1: horizontal profile

        (default = 'vertical')
    image : 2D array-like or None, optional
        The image data which is to be used in the fit. If `None` is specified,
        the function will attempt to grab a new image with the currently
        configured camera at the beamline. (default = None)
    width : int, optional
        Full width of the integration window perpendiculat to the profile
        direction in number of pixels. If None is specified, the width is taken
        to be the full image size. (default = 40)
    trim : int or [int, int], optional
        The number of pixels to trim from the image edge at the start and end
        of the profile. If only one number is specified, the profile is
        trimmed symmetrically. (default = 0)

    Returns
    -------
    fit_result : :class:`lmfit.ModelResult`
        The result of the fit, returned as the corresponding ModelResult
        instance produced by lmfit.

    See also
    --------
    image_fit_profile : Fit a specified profile to image data.

    """

    fit = image_fit_profile(im=image, profile='step',
        background='constant', axis=direction, width=width, trim=trim)

    fit.plot()

    beam_center = fit.params['center'].value
    x = fit.userkws[fit.model.independent_vars[0]]
    cam_center = x[int(len(x)/2.0)]

    y_min, y_max = plt.ylim()
    plt.vlines(cam_center, y_min, y_max, linewidth=1, color='k',
        label=f'det cen: {cam_center:.2f}')
    plt.vlines(beam_center, y_min, y_max, linewidth=2, color='r',
        label=f'edge pos: {beam_center:.2f}')
    plt.legend()

    print("")
    print("Fit result:")
    print("  Detector center position: {} px".format(cam_center))
    print("  Beam center position: {} px".format(beam_center))
    print("")

    return fit



def fit_beam_profile(direction='vertical', image=None, width=40, trim=0):
    """
    Fit a Pseudo-Voigt beam profile on a constant background to the central 40
    pixel wide region of the image along the specified direction.

    Plot the result and display the fitted center position.

    Parameters
    ----------
    direction : str or int
        The axis direction along which the profile is to be computed. Valid
        options are:

        * 'v', 'vertical', or 0: vertical profile
        * 'h', 'horizontal', or 1: horizontal profile

        (default = 'vertical')
    image : 2D array-like or None, optional
        The image data which is to be used in the fit. If `None` is specified,
        the function will attempt to grab a new image with the currently
        configured camera at the beamline. (default = None)
    width : int, optional
        Full width of the integration window perpendiculat to the profile
        direction in number of pixels. If None is specified, the width is taken
        to be the full image size. (default = 40)
    trim : int or [int, int], optional
        The number of pixels to trim from the image edge at the start and end
        of the profile. If only one number is specified, the profile is
        trimmed symmetrically. (default = 0)

    Returns
    -------
    fit_result : :class:`lmfit.ModelResult`
        The result of the fit, returned as the corresponding ModelResult
        instance produced by lmfit.

    See also
    --------
    image_fit_profile : Fit a specified profile to image data.

    """

    fit = image_fit_profile(im=image, profile='pseudo-voigt',
        background='constant', axis=direction, width=width, trim=trim)

    fit.plot()

    beam_center = fit.params['center'].value
    x = fit.userkws[fit.model.independent_vars[0]]
    cam_center = x[int(len(x)/2.0)]

    y_min, y_max = plt.ylim()
    plt.vlines(cam_center, y_min, y_max, linewidth=1, color='k',
        label=f'det cen: {cam_center:.2f}')
    plt.vlines(beam_center, y_min, y_max, linewidth=2, color='r',
        label=f'beam pos: {beam_center:.2f}')
    plt.legend()

    print("")
    print("Fit result:")
    print("  Detector center position: {} px".format(cam_center))
    print("  Beam center position: {} px".format(beam_center))
    print("")

    return fit


def fit_profile(x, y, profile='linear', background=None, plot=False,
                report=False):
    """
    Fit a profile function and optionally a background function to line data.

    Parameters
    ----------
    x : 1D array-like
        The independent variable data of the profile data to be fitted.
    y : 1D array-like
        The dependent variable data of the profile data to be fitted.
    profile : str, optional
        The type of profile to be fitted. The following choices are available:
        'constant', 'linear', 'gaussian', 'lorentzian', 'pseudo-voigt',
        'rectangle', 'step'. (default='linear')
    background : str or None, optional
        The background function to be used in the fit. If `None` or 'none' is
        specified, no background function will be used and only the `profile`
        function is used to fit the data. Valid choices are: 'none',
        'constant', 'linear'. (default = None)
    plot : bool, optional
        If `True`, the function will attempt to plot the result (depends on
        whether matplotlib is available) (default = False)
    report : bool, optional
        If `True`, the fit report will be printed to the command line.
        (default = False)

    Returns
    -------
    fit_result : :class:`lmfit.ModelResult`
        The result of the fit, returned as the corresponding ModelResult
        instance produced by lmfit.

    """

    _profiles = {
        'constant': ConstantModel,
        'linear': LinearModel,
        'gaussian': GaussianModel,
        'lorentzian': LorentzianModel,
        'pseudo-voigt': PseudoVoigtModel,
        'rectangle': RectangleModel,
        'step': StepModel
    }
    _profile_kwargs = {
        'constant': {},
        'linear': {},
        'gaussian': {},
        'lorentzian': {},
        'pseudo-voigt': {},
        'rectangle': {'form': 'erf'},
        'step': {'form': 'erf'}
    }

    _backgrounds = {
        'none': None,
        'constant': ConstantModel,
        'linear': LinearModel,
    }

    if not profile.lower() in _profiles:
        raise ValueError("Not a valid choice for profile: {}\n Valid "
            "choices are {}".format(profile, set(_profiles.keys())))
    if background is None:
        background = 'none'
    if not background.lower() in _backgrounds:
        raise ValueError("Not a valid choice for background: {}\n Valid "
            "choices are {}".format(profile, set(_profiles.keys())))

    x = np.asarray(x)
    y = np.asarray(y)
    if not x.dtype == y.dtype:
        # make sure that both data types are compatible
        # choose float to be on the safe side.
        x = x.astype(float)
        y = y.astype(float)

    prf = _profiles[profile.lower()](**_profile_kwargs[profile.lower()])
    p0 = prf.guess(y - y.min(), x=x)
    model = prf
    if not _backgrounds[background] is None:
        bkgr = _backgrounds[background.lower()]()
        model += bkgr
        p0  += bkgr.guess(np.array([y[0], y[-1]]), x = np.array([x[0], x[-1]]))

    fit_result = model.fit(y, p0, x=x)
    if plot:
        fit_result.plot()
    if report:
        print(fit_result.fit_report())
    return fit_result


def image_fit_profile(im=None, profile='linear', background=None, axis=0,
                      width=None, position=None, trim=0, plot=False,
                      report=False):
    """
    Fit a line profile through (a region of) image data

    Parameters
    ----------
    im : 2D numpy array-like or None, optional
        The image data which is to be used in the fit. If None is specified,
        the function will attempt to grab a new image with the currently
        configured camera at the beamline. (default = None)
    profile : str, optional
        The type of profile to be fitted. The following choices are available:
        'constant', 'linear', 'gaussian', 'lorentzian', 'pseudo-voigt',
        'rectangle', 'step'. (default='linear')
    background : str or None, optional
        The background function to be used in the fit. If None or 'none' is
        specified, no background function will be used and only the `profile`
        function is used to fit the data. Valid choices are: 'none',
        'constant', 'linear'. (default = None)
    axis : str or int, optional
        The axis direction along which the profile is to be computed. Valid
        options are:

        * 'v', 'vertical', or 2: vertical profile
        * 'h', 'horizontal', or 1: horizontal profile

        (default = 0 --> vertical)
    width : int, optional
        Full width of the integration window perpendiculat to the profile
        direction in number of pixels. If None is specified, the width is taken
        to be the full image size. (default = None)
    position : int, optional
        The center position of the profile in pixel coordinates perpendicular
        to the beam direction. If None is specified, the position is set to the
        middle of the image. (default = None)
    trim : int or [int, int], optional
        The number of pixels to trim from the image edge at the start and end
        of the profile. If only one number is specified, the profile is
        trimmed symmetrically. (default = 0)
    plot : bool, optional
        If True, the function will attempt to plot the result (depends on
        whether matplotlib is available) (default = False)
    report : bool, optional
        If True, the fit report will be printed to the command line.
        (default = False)

    Returns
    -------
    fit_result : :class:`lmfit.ModelResult`
        The result of the fit, returned as the corresponding ModelResult
        instance produced by lmfit.

    """

    if im is None:
        # try to grab an image
        bl = get_beamline_config()
        cam = get_camera(bl)
        im = grab_image(camera=cam)

    x,y = get_profile(im, axis=axis, width=width, position=position, trim=trim)
    fit_result = fit_profile(x, y, profile=profile, background=background,
        plot=plot, report=report)

    return fit_result
