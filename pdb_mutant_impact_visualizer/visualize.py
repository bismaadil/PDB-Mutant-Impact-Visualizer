"""3D C-alpha view of an input structure, colored by deposited B-factor."""

import plotly.graph_objects as go
from Bio.PDB.Structure import Structure

from .local_structure import THREE_TO_ONE, analyze_local_structure, get_residue


def build_structure_figure(
    structure: Structure, chain_id: str, position: int, cutoff: float = 10.0,
    insertion_code: str = " ",
) -> go.Figure:
    target = get_residue(structure, chain_id, position, insertion_code)
    local = analyze_local_structure(structure, chain_id, position, cutoff, insertion_code)
    neighbor_keys = {
        (item["chain"], item["author_residue_number"], item["insertion_code"])
        for item in local["neighbors"]
    }
    ca_atoms = [r["CA"] for r in structure[0][chain_id] if "CA" in r and THREE_TO_ONE.get(r.resname)]
    xs, ys, zs, colors, labels, sizes = [], [], [], [], [], []
    for atom in ca_atoms:
        residue = atom.get_parent()
        xs.append(float(atom.coord[0])); ys.append(float(atom.coord[1])); zs.append(float(atom.coord[2]))
        colors.append(float(atom.bfactor))
        labels.append(f"{chain_id}:{residue.id[1]}{residue.id[2].strip()} {THREE_TO_ONE[residue.resname]} | B={atom.bfactor:.1f} Å²")
        is_neighbor = (chain_id, residue.id[1], residue.id[2].strip()) in neighbor_keys
        label = "MUTATION SITE" if residue is target else "NEIGHBOR" if is_neighbor else ""
        if label:
            labels[-1] = f"{label} | " + labels[-1]
        sizes.append(13 if residue is target else 8 if is_neighbor else 4)
    figure = go.Figure(go.Scatter3d(
        x=xs, y=ys, z=zs, mode="lines", line={"color": "rgba(120,120,120,0.45)", "width": 3},
        hoverinfo="skip", showlegend=False,
    ))
    figure.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs, mode="markers", text=labels, hoverinfo="text",
        marker={"size": sizes, "color": colors, "colorscale": "Viridis", "colorbar": {"title": "B-factor (Å²)"}, "opacity": 0.9}, showlegend=False,
    ))
    figure.update_layout(
        title=f"Chain {chain_id} Cα atoms — mutation site {position}",
        scene={"xaxis_title": "X (Å)", "yaxis_title": "Y (Å)", "zaxis_title": "Z (Å)", "aspectmode": "data"},
        margin={"l": 0, "r": 0, "t": 45, "b": 0},
    )
    return figure
