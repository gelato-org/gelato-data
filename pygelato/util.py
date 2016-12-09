# coding=utf8
from __future__ import unicode_literals

from clldutils.path import Path

import pygelato

REPOS_PATH = Path(pygelato.__file__).parent.parent


def data_path(*comps, **kw):
    return kw.get('repos', REPOS_PATH).joinpath('datasets', *comps)
