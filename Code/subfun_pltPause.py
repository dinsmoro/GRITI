#8b503459569d1505560191b614e39403357c9bce68a28ae15a01700b47aae0c888a8a51f3779cce626892184be3565f5207cc09f4cef1cf8c97613e0fbf87963
import time
from matplotlib import _pylab_helpers
def pause(interval):
    """
    Run the GUI event loop for *interval* seconds.

    If there is an active figure, it will be updated and displayed before the
    pause, and the GUI event loop (if any) will run during the pause.

    This can be used for crude animation.  For more complex animation use
    :mod:`matplotlib.animation`.

    If there is no active figure, sleep for *interval* seconds instead.

    See Also
    --------
    matplotlib.animation : Proper animations
    show : Show all figures and optional block until all figures are closed.
    """
    manager = _pylab_helpers.Gcf.get_active()
    if manager is not None:
        canvas = manager.canvas
        if canvas.figure.stale:
            canvas.draw_idle()
        canvas.start_event_loop(interval)
    else:
        time.sleep(interval)


