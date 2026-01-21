from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from zoneinfo import ZoneInfo


@dataclass(frozen=True)
class RunTime:
    date: str  # YYYY-MM-DD
    time: str  # HHMMSS


def now_jst() -> RunTime:
    dt = datetime.now(ZoneInfo("Asia/Tokyo"))
    return RunTime(date=dt.strftime("%Y-%m-%d"), time=dt.strftime("%H%M%S"))
