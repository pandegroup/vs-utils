"""
Featurize molecules and save features to disk. Featurizers are exposed as
subcommands, with __init__ arguments as subcommand arguments.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import cPickle
import inspect

from pande_gas.features import get_featurizers
from pande_gas.utils import h5
from pande_gas.utils import rdkit_utils as rd


class HelpFormatter(argparse.RawTextHelpFormatter):
    """
    Argparse help formatter with better indenting.

    Parameters
    ----------
    WRITEME
    """
    def __init__(self, prog, indent_increment=2, max_help_position=8,
                 width=None):
        super(HelpFormatter, self).__init__(prog, indent_increment,
                                            max_help_position, width)


def main():
    mols = rd.read(args.input)
    names = [mol.GetProp('_Name') for mol in mols]
    featurizer = args.klass(**vars(args.featurizer_kwargs))
    features = featurizer.featurize(mols)
    data = {'features': features, 'names': names}
    if args.labels is not None:
        with open(args.labels) as f:
            labels = cPickle.load(f)
        assert len(labels) == len(mols)
        data['y'] = labels
    h5.dump(data, args.output, attrs=vars(args.featurizer_kwargs))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument('input',
                        help='Input molecules.')
    parser.add_argument('-l', '--labels',
                        help='Molecule labels.')
    parser.add_argument('output',
                        help='Output (HDF5) filename.')

    # featurizer subcommands
    featurizers = get_featurizers()
    subparsers = parser.add_subparsers(title='featurizers')
    for name, klass in featurizers.items():
        command = subparsers.add_parser(name, help=klass.__doc__,
                                        formatter_class=HelpFormatter,
                                        epilog=klass.__doc__)
        command.set_defaults(klass=klass)
        try:
            args, _, _, defaults = inspect.getargspec(klass.__init__)
        except TypeError:
            args = []
        for i, arg in enumerate(args):
            if i == 0 and arg == 'self':
                continue
            kwargs = {}
            try:
                kwargs['default'] = defaults[i-len(args)]
                kwargs['type'] = type(kwargs['default'])
            except IndexError:
                kwargs['required'] = True
            if 'type' in kwargs and kwargs['type'] == bool:
                if kwargs['default']:
                    command.add_argument('--no-{}'.format(arg), dest=arg,
                                         action='store_false')
                else:
                    command.add_argument('--{}'.format(arg),
                                         action='store_true')
            else:
                command.add_argument('--{}'.format(arg), **kwargs)
    args = argparse.Namespace()
    args.featurizer_kwargs = parser.parse_args()
    for arg in ['input', 'output', 'klass', 'labels']:
        setattr(args, arg, getattr(args.featurizer_kwargs, arg))
        delattr(args.featurizer_kwargs, arg)
    main()
