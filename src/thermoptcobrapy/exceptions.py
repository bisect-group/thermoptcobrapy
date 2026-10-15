"""The errors ThermOpt raises."""


class ThermOptError(Exception):
    """Base class for every error this package raises."""


class ThermOptSolverError(ThermOptError):
    """A solver did not reach optimality when it had to.

    Parameters
    ----------
    message : str
        What went wrong.
    status : str, optional
        The solver status that was reached instead (default None).

    """

    def __init__(self, message: str, status: str = None) -> None:
        super().__init__(message)
        self.status = status


class ThermOptInfeasibleCoreError(ThermOptError):
    """No thermodynamically feasible flux lets the core carry flux at once."""


class ThermOptTimeLimitError(ThermOptError):
    """A mixed-integer problem hit its time limit with nothing to show."""
