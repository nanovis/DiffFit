import numpy as np


def shift_difference(shift1, shift2):
    """
    Calculate the Euclidean distance between two shifts.
    """
    return np.linalg.norm(shift1 - shift2)


def quaternion_angle_distance(q1, q2):
    """
    Calculate the angular distance in degrees between two quaternions.
    """
    q1 = q1 / np.linalg.norm(q1)
    q2 = q2 / np.linalg.norm(q2)
    dot = np.dot(q1, q2)
    dot = np.clip(dot, -1.0, 1.0)
    angle_rad = 2 * np.arccos(abs(dot))
    angle_deg = np.degrees(angle_rad)
    return angle_deg


def test_rl():
    print("Test RL")


def cluster_sqd():
    return 0


def sort_sqd_cluster():
    return 0