from ase.io.trajectory import Trajectory
from ase.io import read, write

def split_trajectory(input_file, output_prefix, frames_per_file=100):
    """
    Splits a large ASE trajectory file into smaller files.

    Args:
        input_file (str): Path to the input trajectory file.
        output_prefix (str): Prefix for the output filenames (e.g., "traj_").
        frames_per_file (int): Number of frames to include in each output file.
    """
    traj = Trajectory(input_file)
    num_frames = len(traj)
    frames_per_file = int(num_frames/2)
    output_index = 1

    for i in range(0, num_frames, frames_per_file):
        output_file = f"{output_prefix}{output_index:04d}.traj"
        with Trajectory(output_file, 'w') as output_traj:
            for j in range(i, min(i + frames_per_file, num_frames)):
                atoms = traj[j]
                output_traj.write(atoms)
        output_index += 1
    traj.close()


# Example usage:
input_trajectory_file = "md.traj"
output_file_prefix = "split_md"
split_trajectory(input_trajectory_file, output_file_prefix, frames_per_file=1)
