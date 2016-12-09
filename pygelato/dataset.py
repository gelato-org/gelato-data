# coding: utf8
from __future__ import unicode_literals, print_function, division
import logging

from clldutils import jsonlib
from clldutils.path import Path, import_module

import pygelato
from pygelato.util import data_path


logging.basicConfig(level=logging.INFO)


class Dataset(object):
    """
    A gelato dataset.

    This object provides access to a dataset's
    - python module as attribute `commands`
    - language list as attribute `languages`
    """
    def __init__(self, path, glottolog_repos=None):
        """
        A dataset is initialzed by passing its directory path.
        """
        path = Path(path)
        self.id = path.name
        self.log = logging.getLogger(pygelato.__name__)
        self.dir = path

        # raw data
        self.raw = self.dir.joinpath('raw', 'data')
        if not self.raw.exists():
            self.raw.mkdir()

        # processed data
        self.processed = self.dir.joinpath('processed')
        if not self.processed.exists():
            self.processed.mkdir()

        self.commands = import_module(self.dir)
        self.md = jsonlib.load(self.dir.joinpath('metadata.json'))

    @classmethod
    def from_name(cls, name):
        """
        Factory method, looking up a dataset by name in the default data directory.
        """
        return cls(data_path(name))

    def _run_command(self, name, *args, **kw):
        """
        Call a callable defined in the dataset's python module, if available.

        :param name: Name of the callable.
        :param args: Positional arguments are passed to the callable.
        :param kw: Keyword arguments are passed to the callable.
        :return:
        """
        if not hasattr(self.commands, name):
            self.log.warn('command "%s" not available for dataset %s' % (name, self.id))
        else:
            getattr(self.commands, name)(self, *args, **kw)

    def download(self, **kw):
        self._run_command('download', **kw)
