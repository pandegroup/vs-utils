"""
Ligand and receptor preparation for AutoDock/Vina.
"""
import argparse
import os
import subprocess
import tempfile


def parse_args(input_args=None):
  """
  Parse command-line arguments.

  Parameters
  ----------
  input_args : list, optional
    Input arguments. If not provided, defaults to sys.argv[1:].
  """
  parser = argparse.ArgumentParser()
  parser.add_argument('input',
                      help='Input filename.')
  parser.add_argument('-r', '--receptor', action='store_true',
                      help='Treat input as a receptor.')
  parser.add_argument('-s', '--split', action='store_true',
                      help='Split output files into one file per molecule.')
  return parser.parse_args(input_args)


def main(input_filename, receptor=False, split=False):
  output_filename = os.path.basename(input_filename)
  if output_filename.endswith('.gz'):
    output_filename = '.'.join(output_filename.split('.')[:-2])
  else:
    output_filename = '.'.join(output_filename.split('.')[:-1])
  output_filename += '.pdbqt'
  if receptor:
    assert input_filename.endswith('.pdb'), 'Input file must be a PDB.'

    # process PDB with pdbfixer
    f, fixed_filename = tempfile.mkstemp(suffix='-{}'.format(input_filename))
    os.close(f)
    try:
      pdbfixer_args = ['pdbfixer', input_filename, '--keep-heterogens=none',
                       '--ph=7.4', '--add-atoms=hydrogen',
                       '--output={}'.format(fixed_filename)]
      subprocess.check_call(pdbfixer_args)
      obabel_args = ['obabel', fixed_filename, '-O', output_filename, '-xr']
      subprocess.check_call(obabel_args)
    finally:
      os.remove(fixed_filename)
  else:
    obabel_args = ['obabel', input_filename, '-O', output_filename]
    if split:
      obabel_args += ['-m']
    subprocess.check_call(obabel_args)
  print '{} -> {}'.format(input_filename, output_filename)

if __name__ == '__main__':
  args = parse_args()
  main(args.input, args.receptor, args.split)
