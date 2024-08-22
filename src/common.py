import torch
import matplotlib.pyplot as plt
import numpy as np

def generate_grid_coordinates(start: float=-0.03, end: float=0.03, steps: int=3):
    # Example grid size for each dimension
    # To make the grid size controllable, adjust steps for x, y, and z
    steps_x = steps
    steps_y = steps
    steps_z = steps

    # Create linear spaces for each dimension
    x = torch.linspace(start, end, steps=steps_x)
    y = torch.linspace(start, end, steps=steps_y)
    z = torch.linspace(start, end, steps=steps_z)

    # Create a meshgrid for the coordinates
    X, Y, Z = torch.meshgrid(x, y, z, indexing='ij')

    # Combine the coordinates into a single tensor of shape (steps_x, steps_y, steps_z, 3)
    grid_coordinates = torch.stack([X, Y, Z], dim=-1)
    return grid_coordinates


# from DiFit.common import histogram_tensor
def histogram_tensor(tensor):
    tensor_flattened = tensor.flatten().cpu().numpy()
    plt.figure(figsize=(10, 6))
    plt.hist(tensor_flattened, bins=50, alpha=0.75, color='blue')
    plt.title('Histogram of Tensor Values')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()
    plt.close()


def generate_random_quaternions(n):
    """
    Generate n random quaternion vectors that evenly sample directions.

    Parameters:
    - n: The number of quaternion vectors to generate.

    Returns:
    - quaternions: An array of shape (n, 4) containing n quaternions.
    """
    # Random rotation axes
    u1 = np.random.rand(n)
    u2 = np.random.rand(n)
    u3 = np.random.rand(n)

    # Convert random numbers to quaternion parameters
    q0 = np.sqrt(1 - u1) * np.sin(2 * np.pi * u2)
    q1 = np.sqrt(1 - u1) * np.cos(2 * np.pi * u2)
    q2 = np.sqrt(u1) * np.sin(2 * np.pi * u3)
    q3 = np.sqrt(u1) * np.cos(2 * np.pi * u3)

    # Combine into quaternion array
    quaternions = np.vstack((q0, q1, q2, q3)).T

    return quaternions


def quaternion_angle_distance(q1, q2):
    # Function to calculate the angular distance in degrees between two quaternions

    # Ensure quaternions are normalized
    q1 = q1 / np.linalg.norm(q1)
    q2 = q2 / np.linalg.norm(q2)

    # Dot product between two quaternions
    dot = np.dot(q1, q2)

    # Clamp dot product to ensure it's within the valid range for arccos
    dot = np.clip(dot, -1.0, 1.0)

    # Calculate the angle in radians and then convert to degrees
    angle_rad = 2 * np.arccos(abs(dot))
    angle_deg = np.degrees(angle_rad)

    return angle_deg