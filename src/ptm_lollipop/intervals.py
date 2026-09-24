from __future__ import annotations

import re


def parse_ranges(text: str, length: int | None = None) -> tuple[tuple[int, int], ...]:
    value = (text or "").strip().lower()
    if not value or value in {"no", "0", "(0)", "none", "na", "n/a"}:
        return ()
    if value == "yes":
        return ()

    ranges: list[tuple[int, int]] = []
    for start_s, end_s in re.findall(r"(\d+)\s*[-:]\s*(\d+)", value):
        start, end = int(start_s), int(end_s)
        if length is not None:
            start, end = max(1, start), min(length, end)
        if end >= start:
            ranges.append((start, end))
    return tuple(ranges)


def merge_intervals(intervals: list[tuple[int, int]] | tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    if not intervals:
        return ()
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return tuple(merged)


def fmt_ranges(intervals: tuple[tuple[int, int], ...] | list[tuple[int, int]]) -> str:
    return "; ".join(f"{start}-{end}" for start, end in intervals) if intervals else "none"


def interval_overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def compare_disorder(user_text: str, sgd_intervals: tuple[tuple[int, int], ...], length: int) -> tuple[str, str]:
    value = (user_text or "").strip().lower()
    user = parse_ranges(user_text, length)
    if value == "yes":
        return "user says yes, no coordinates", "manual review"

    user_aa = sum(end - start + 1 for start, end in user)
    sgd_aa = sum(end - start + 1 for start, end in sgd_intervals)
    overlap = sum(interval_overlap(a, b) for a in user for b in sgd_intervals)

    if not user and not sgd_intervals:
        return "matches no disorder", "ok"
    if not user and sgd_intervals:
        return f"user none/blank; SGD {fmt_ranges(sgd_intervals)}", "SGD differs"
    if user and not sgd_intervals:
        return f"user {fmt_ranges(user)}; SGD none", "SGD differs"
    if user == sgd_intervals:
        return "exact match", "ok"

    recall = overlap / user_aa if user_aa else 1.0
    precision = overlap / sgd_aa if sgd_aa else 1.0
    if recall >= 0.85 and precision >= 0.85:
        return f"similar; user {fmt_ranges(user)} vs SGD {fmt_ranges(sgd_intervals)}", "near match"
    return f"differs; user {fmt_ranges(user)} vs SGD {fmt_ranges(sgd_intervals)}", "SGD differs"

