from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable, List, TypeVar

T = TypeVar("T")
R = TypeVar("R")


def run_threaded(
    fn: Callable[[T], R],
    items: Iterable[T],
    threads: int,
) -> List[R]:
    """
    Run fn(item) for each item in items, possibly in parallel.
    """
    if threads <= 1:
        return [fn(x) for x in items]

    with ThreadPoolExecutor(max_workers=threads) as ex:
        return list(ex.map(fn, items))
