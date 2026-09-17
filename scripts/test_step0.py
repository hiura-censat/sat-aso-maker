import gzip, tempfile, unittest, csv
from pathlib import Path
import step0_preprocess as p
class Regions(unittest.TestCase):
    def test_overlap_partition(self):
        a=p.merge([(2,10),(8,20),(20,22)]); b=[(18,30)]
        overlap=p.intersect(a,b)
        self.assertEqual(a,[(2,22)])
        self.assertEqual(overlap,[(18,22)])
        self.assertEqual(p.subtract(a,overlap),[(2,18)])
        self.assertEqual(p.subtract([(0,40)],p.merge(a+b)),[(0,2),(30,40)])
    def test_windows_ambiguity_case_and_boundaries(self):
        with tempfile.TemporaryDirectory() as d:
            old=p.OUT; p.OUT=Path(d)
            try:
                fa=Path(d)/'input.gz'
                seq=b'A'*16+b'N'+b'c'*17+b'R'+b'T'*15
                with gzip.open(str(fa),'wb') as f: f.write(b'>chr1\n'+seq+b'\n')
                p.stats(('test',str(fa),{'chr1':len(seq)},{'whole':{'chr1':[(0,len(seq))]},'split':{'chr1':[(0,8),(8,16)]}}))
                with (Path(d)/'qc/test.stats.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
                self.assertEqual(int(rows[0]['valid_16mer_starts']),3)
                self.assertEqual(int(rows[0]['non_ACGT_bp']),2)
                self.assertEqual(int(rows[1]['valid_16mer_starts']),0)
            finally: p.OUT=old
if __name__=='__main__': unittest.main()
