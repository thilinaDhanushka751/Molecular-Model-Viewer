from flask import Flask, request, jsonify, render_template
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/generate_3d', methods=['POST'])
def generate_3d():
    data = request.get_json()
    formula = data.get('formula')

    try:
        mol = Chem.MolFromSmiles(formula)
        if mol is None:
            return jsonify({'error': 'Invalid chemical notation'}), 400

        # Add hydrogens and embed molecule
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            return jsonify({'error': 'Failed to generate 3D structure'}), 400

        # Optimize geometry
        if AllChem.UFFOptimizeMolecule(mol) != 0:
            return jsonify({'error': 'Failed to optimize molecule geometry'}), 400

        # Extract 3D coordinates
        atoms = []
        for atom in mol.GetAtoms():
            atoms.append({
                'symbol': atom.GetSymbol(),
                'index': atom.GetIdx()
            })

        bonds = []
        for bond in mol.GetBonds():
            # Map bond types to order (1 = single, 2 = double, 3 = triple)
            bond_order = {
                Chem.BondType.SINGLE: 1,
                Chem.BondType.DOUBLE: 2,
                Chem.BondType.TRIPLE: 3,
                Chem.BondType.AROMATIC: 1.5  # Aromatic bonds as 1.5 (optional, for visualization)
            }.get(bond.GetBondType(), 1)  # Default to single if not in the map

            bonds.append({
                'start': bond.GetBeginAtomIdx(),
                'end': bond.GetEndAtomIdx(),
                'order': bond_order
            })

        positions = []
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            positions.append([pos.x, pos.y, pos.z])

        return jsonify({
            'atoms': atoms,
            'bonds': bonds,
            'positions': positions
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)
