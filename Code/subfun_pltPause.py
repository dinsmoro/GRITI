#8777e25fc1e33a8b23dd3b8ae0cae99a65d45f20f5cd8e174fcd02a89292f463bab0e1cc24c4d65c2dbb3d13f03e75791b89b62baf7b8c18938a4c98a1badcbd
import time
from matplotlib import _pylab_helpers
def pause(interval: float) -> None:
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


