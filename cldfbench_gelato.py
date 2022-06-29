import pathlib

from cldfbench import Dataset as BaseDataset


class Dataset(BaseDataset):
    dir = pathlib.Path(__file__).parent
    id = "gelato"

    def cldf_specs(self):  # A dataset must declare all CLDF sets it creates.
        return super().cldf_specs()

    def cmd_download(self, args):
        """
        Download files to the raw/ directory. You can use helpers methods of `self.raw_dir`, e.g.

        >>> self.raw_dir.download(url, fname)
        """
        pass

    def cmd_makecldf(self, args):
        for d in self.dir.joinpath('datasets').iterdir():
            if d.is_dir() and d.stem != 'Pemberton_AutosomalSTR':
                for s in d.read_csv('samples.csv', dicts=True):
                    print(list(s.keys()))
                    break
