from pathlib import Path
from monty.serialization import loadfn as _loadfn
from ase.io import Trajectory as AseTrajectory
from pymatgen.io.ase import AseAtomsAdaptor
from tempfile import NamedTemporaryFile
from ase.io import write
from sys import argv
from joblib import Memory

memory = Memory(".cache", verbose=0)

def loadfn(path: str | Path):
    """
    Load a file using monty serialization.

    Args:
        path (str or Path): The path to the file to load.

    Returns:
        The loaded object.
    """
    return memory.cache(_loadfn)(path)

def to_ase(
    traj_pmg,
    ase_traj_file: str | Path | None = None,
    oxi_state_map: dict[str, int] | None = None
) -> AseTrajectory:
    """
    Convert a pymatgen .Trajectory to an ASE .Trajectory.

    Args:
        trajectory (pymatgen .Trajectory) : trajectory to convert
        property_map (dict[str,str]) : A mapping between ASE calculator properties and
            pymatgen .Trajectory `frame_properties` keys. Ex.:
                property_map = {"energy": "e_0_energy"}
            would map `e_0_energy` in the pymatgen .Trajectory `frame_properties`
            to ASE's `get_potential_energy` function.
            See `ase.calculators.calculator.all_properties` for a list of acceptable calculator properties.
        ase_traj_file (str, Path, or None (default) ) :  If not None, the name of
            the file to write the ASE trajectory to.

    Returns:
        ase .Trajectory
    """

    adaptor = AseAtomsAdaptor()
    frame_props = traj_pmg.frame_properties or [{} for _ in range(len(traj_pmg))]
    temp_file = None
    if ase_traj_file is None:
        temp_file = NamedTemporaryFile(delete=False)  # noqa: SIM115
        ase_traj_file = temp_file.name

    for idx, structure in enumerate(traj_pmg):
        atoms = adaptor.get_atoms(structure, msonable=False, velocities=structure.site_properties.get("velocities"))



        atoms.info["REF_energy"] = frame_props[idx]['e_0_energy']
        atoms.info["REF_stress"] = frame_props[idx]['stress']
        atoms.arrays["REF_forces"] = frame_props[idx]['forces']
        if oxi_state_map is not None:
            atoms.arrays['REF_formal_charge'] = [oxi_state_map[symbol] for symbol in atoms.get_chemical_symbols()]


        with AseTrajectory(ase_traj_file, "a" if idx > 0 else "w", atoms=atoms) as _traj_file:
            _traj_file.write()

    ase_traj = AseTrajectory(ase_traj_file, "r")
    if temp_file is not None:
        temp_file.close()

    return ase_traj


DATA_DIR_PATH = Path("/global/cfs/projectdirs/matgen/virkaran/MD_trajectories/BTO_full_prod_traj/")
OXI_STATE_MAP = {"O": -2, "Ti": 4, "Ba": 2}
_, RANK = argv
RANK = int(RANK)


paths = DATA_DIR_PATH.iterdir()
paths = sorted(paths)
path = paths[RANK]

traj_name = path.stem
traj = loadfn(path)
ase_traj = to_ase(traj, ase_traj_file=Path("data_traj")/f"{traj_name}.traj", oxi_state_map=OXI_STATE_MAP)
write(Path("data_xyz")/f"{traj_name}.xyz", ase_traj, format="xyz")

