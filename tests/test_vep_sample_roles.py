import ast
from pathlib import Path
import unittest

# Test the identity function from the production source without optional annotation databases.
source = Path(__file__).resolve().parents[1] / 'scripts' / 'vep_report.py'
node = next(n for n in ast.parse(source.read_text(encoding='utf-8')).body if isinstance(n, ast.FunctionDef) and n.name == 'align_sample_columns')
namespace = {}
exec(compile(ast.Module(body=[node], type_ignores=[]), str(source), 'exec'), namespace)
align = namespace['align_sample_columns']

class SampleRolesTest(unittest.TestCase):
    def test_vcf_order_is_not_report_order(self):
        original = {'sample_names':['father','mother','proband'], 'samples':[{'gt':'0/0','ad':'53,0'}, {'gt':'0/1','ad':'58,59'}, {'gt':'1|0','ad':'63,68'}]}
        result = align([original], ['proband','father','mother'])[0]
        self.assertEqual([s['ad'] for s in result['samples']], ['63,68','53,0','58,59'])
        self.assertEqual(original['samples'][0]['ad'], '53,0')
        self.assertEqual(align([result], ['proband','father','mother']), [result])
    def test_unknown_and_duplicate_identities_fail(self):
        for names in [['father','father','proband'], [], ['sample1','sample2','sample3']]:
            with self.assertRaises(ValueError):
                align([{'sample_names':names,'samples':[{}, {}, {}]}], ['proband','father','mother'])
    def test_single_sample_alias_remains_supported(self):
        self.assertEqual(align([{'sample_names':['external'],'samples':[{'gt':'0/1'}]}],['proband'])[0]['samples'],[{'gt':'0/1'}])

if __name__ == '__main__': unittest.main()
