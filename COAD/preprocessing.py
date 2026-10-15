import pandas as pd
import numpy as np
import warnings
import os
warnings.filterwarnings('ignore')

class Preprocessing:
    def __init__(self, miRNACancer, miRNACancerlength, miRNANormal, miRNANormallength, geneCancer, geneNormal, tumor):
        self.miRNACancer = miRNACancer.copy()
        self.miRNACancerlength = miRNACancerlength.copy()
        self.miRNANormal = miRNANormal.copy()
        self.miRNANormallength = miRNANormallength.copy()
        self.geneCancer = geneCancer.copy()
        self.geneNormal = geneNormal.copy()
        self.tumor = tumor

        self.fix_matrices()
        self.geneCancer, self.geneNormal = self.filter_genes()
        self.filter_noise()
        self.common_patients()
        self.Normal, self.Cancer, self.Indices = self.matrices()

    def fix_matrices(self):
        for df in [self.miRNACancer, self.miRNANormal, self.miRNANormallength, self.miRNACancerlength]:
            if isinstance(df.index, pd.MultiIndex):
                df.index = df.index.droplevel(0)
            else:
                df.index.name = 'gene_symbol'

        for df in [self.miRNANormal, self.miRNACancer]:
            if isinstance(df.columns, pd.MultiIndex):
                df.columns = df.columns.droplevel(0)
        self.miRNANormal.columns.name = None
        self.miRNACancer.columns.name = None

        if isinstance(self.geneNormal.columns, pd.MultiIndex):
            self.geneNormal.columns = self.geneNormal.columns.droplevel(0)
        if isinstance(self.geneCancer.columns, pd.MultiIndex):
            self.geneCancer.columns = self.geneCancer.columns.droplevel(0)
        self.geneNormal.columns.name = None
        self.geneCancer.columns.name = None

        print('geneNorm ' + str(self.geneNormal.shape))
        print('geneCanc ' + str(self.geneCancer.shape))

        l = len(self.miRNANormal)

        if hasattr(self.miRNANormallength, 'shape') and self.miRNANormallength.shape[1] >= 2 and self.miRNANormallength.shape != self.miRNANormal.shape:
            for i in range(0, l - 1):
                self.miRNANormal.values[i] = (self.miRNANormal.values[i] * (10**3)) / (2 * abs(self.miRNANormallength.values[i, 1] - self.miRNANormallength.values[i, 0]))
                self.miRNACancer.values[i] = (self.miRNACancer.values[i] * (10**3)) / (2 * abs(self.miRNACancerlength.values[i, 1] - self.miRNACancerlength.values[i, 0]))

        self.miRNANormal.index.name = 'ensembl'
        self.miRNACancer.index.name = 'ensembl'

    def filter_genes(self):
        filename = 'matrix_GTypeGName.xls'
        if not os.path.exists(filename):
            raise FileNotFoundError(f"❌ The required annotation file '{filename}' is missing.")

        matrix_GTypeGName = pd.read_excel(filename)

        # 1. Recover the patient column if it was accidentally pushed into the index
        if isinstance(self.geneCancer.index, pd.MultiIndex):
            self.geneCancer = self.geneCancer.reset_index(level=1)
        if isinstance(self.geneNormal.index, pd.MultiIndex):
            self.geneNormal = self.geneNormal.reset_index(level=1)

        # 2. Extract clean Ensembl IDs
        self.geneCancer['ensembl'] = self.geneCancer.index.astype(str).str.split('.').str[0]
        self.geneNormal['ensembl'] = self.geneNormal.index.astype(str).str.split('.').str[0]

        # 3. Use mygene to translate Ensembl IDs to Gene Symbols dynamically
        print("🔍 Automatically translating Ensembl IDs to Gene Names... (Takes ~10 seconds)")
        try:
            import mygene
        except ImportError:
            raise ImportError("🚨 The 'mygene' library is missing! Run '!pip install mygene' in your notebook first.")

        mg = mygene.MyGeneInfo()
        ensembl_list = self.geneCancer['ensembl'].dropna().unique().tolist()

        # Fetch the mapping dictionary from the internet
        mapping = mg.querymany(ensembl_list, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)

        # Format the translation table
        if 'symbol' in mapping.columns:
            mapping = mapping.reset_index()[['query', 'symbol']].rename(columns={'query': 'ensembl', 'symbol': 'gene_symbol'})
            mapping = mapping.dropna().drop_duplicates(subset='ensembl')
        else:
            mapping = pd.DataFrame(columns=['ensembl', 'gene_symbol'])

        # Attach the translated text symbols to your gene matrices
        self.geneCancer = pd.merge(self.geneCancer, mapping, on='ensembl', how='left')
        self.geneNormal = pd.merge(self.geneNormal, mapping, on='ensembl', how='left')

        # 4. NOW we can safely merge with the Excel file! Both sides are text now.
        GeneCancerTypes = pd.merge(self.geneCancer, matrix_GTypeGName, left_on="gene_symbol", right_on="Gene name")
        GeneNormalTypes = pd.merge(self.geneNormal, matrix_GTypeGName, left_on="gene_symbol", right_on="Gene name")

        GeneCancerTypes.index = GeneCancerTypes['ensembl']
        GeneNormalTypes.index = GeneNormalTypes['ensembl']
        GeneCancerTypes = GeneCancerTypes.drop_duplicates(subset=GeneCancerTypes.columns.difference(['ensembl', 'gene_symbol', 'Gene name', 'Gene type']))
        GeneNormalTypes = GeneNormalTypes.drop_duplicates(subset=GeneNormalTypes.columns.difference(['ensembl', 'gene_symbol', 'Gene name', 'Gene type']))

        # 5. Apply the biological filtering rules
        types_to_exclude = [
            'TR_J_gene', 'polymorphic_pseudogene', 'unitary_pseudogene', 'IG_D_gene', 'scaRNA',
            'TR_V_pseudogene', 'sense_overlapping', 'rRNA', 'IG_J_gene', 'IG_C_gene',
            'bidirectional_promoter_lncRNA', 'pseudogene', 'IG_C_pseudogene', '3prime_overlapping_ncRNA',
            'TR_C_gene', 'IG_J_pseudogene', 'ribozyme', 'TR_J_pseudogene', 'TR_D_gene', 'Mt_rRNA',
            'antisense', 'misc_RNA', 'rRNA_pseudogene', 'transcribed_unprocessed_pseudogene',
            'processed_transcript', 'transcribed_processed_pseudogene', 'IG_V_pseudogene',
            'sense_intronic', 'IG_V_gene', 'TR_V_gene', 'TEC', 'transcribed_unitary_pseudogene',
            'vaultRNA', 'macro_lncRNA', 'scRNA', 'non_coding', 'miRNA', 'snRNA', 'snoRNA', 'Mt_tRNA'
        ]

        for g_type in types_to_exclude:
            GeneCancerTypes = GeneCancerTypes[GeneCancerTypes['Gene type'] != g_type]
            GeneNormalTypes = GeneNormalTypes[GeneNormalTypes['Gene type'] != g_type]

        GeneCancerLong = GeneCancerTypes.drop(['Gene type', 'Gene name', 'ensembl', 'gene_symbol'], axis=1, errors='ignore')
        GeneNormalLong = GeneNormalTypes.drop(['Gene type', 'Gene name', 'ensembl', 'gene_symbol'], axis=1, errors='ignore')

        GeneCancerLong.index.name = 'ensembl'
        GeneNormalLong.index.name = 'ensembl'

        indexDuplicated_cancer = GeneCancerLong[GeneCancerLong.index.duplicated()].index
        indexDuplicated_normal = GeneNormalLong[GeneNormalLong.index.duplicated()].index
        indexDuplicated = indexDuplicated_cancer.intersection(indexDuplicated_normal)

        noDuplCancer = GeneCancerLong.drop(indexDuplicated, errors='ignore')
        noDuplNormal = GeneNormalLong.drop(indexDuplicated, errors='ignore')
        DuplCancer = GeneCancerLong.loc[indexDuplicated].drop_duplicates() if not indexDuplicated.empty else pd.DataFrame()
        DuplNormal = GeneNormalLong.loc[indexDuplicated].drop_duplicates() if not indexDuplicated.empty else pd.DataFrame()

        geneCancer_noDupl = pd.concat([noDuplCancer, DuplCancer], axis=0)
        geneNormal_noDupl = pd.concat([noDuplNormal, DuplNormal], axis=0)

        diff = geneCancer_noDupl.index.difference(geneNormal_noDupl.index)
        geneCancer_noDupl = geneCancer_noDupl.drop(diff, axis=0)
        diff_norm = geneNormal_noDupl.index.difference(geneCancer_noDupl.index)
        geneNormal_noDupl = geneNormal_noDupl.drop(diff_norm, axis=0)

        print('normal ' + str(geneNormal_noDupl.shape))
        print('cancer ' + str(geneCancer_noDupl.shape))
        return geneCancer_noDupl, geneNormal_noDupl

    def filter_noise(self):
        miRNANormal_mean = self.miRNANormal.mean(axis=1)
        miRNACancer_mean = self.miRNACancer.mean(axis=1)
        geneNormal_mean = self.geneNormal.mean(axis=1)
        geneCancer_mean = self.geneCancer.mean(axis=1)

        miRNANormalIndex_null = miRNANormal_mean[miRNANormal_mean.values == 0].index
        miRNACancerIndex_null = miRNACancer_mean[miRNACancer_mean.values == 0].index
        geneNormalIndex_null = geneNormal_mean[geneNormal_mean.values == 0].index
        geneCancerIndex_null = geneCancer_mean[geneCancer_mean.values == 0].index

        miRNA_null = miRNANormalIndex_null.intersection(miRNACancerIndex_null)
        gene_null = geneNormalIndex_null.intersection(geneCancerIndex_null)

        self.miRNANormal = self.miRNANormal.drop(miRNA_null, axis=0, errors='ignore')
        self.miRNACancer = self.miRNACancer.drop(miRNA_null, axis=0, errors='ignore')
        self.geneNormal = self.geneNormal.drop(gene_null, axis=0, errors='ignore')
        self.geneCancer = self.geneCancer.drop(gene_null, axis=0, errors='ignore')

        geneNormal_mean = self.geneNormal.mean(axis=1)
        geneCancer_mean = self.geneCancer.mean(axis=1)

        geneNormal_log2 = np.log2(self.geneNormal.replace(0, np.nan))
        geneCancer_log2 = np.log2(self.geneCancer.replace(0, np.nan))

        geneNormal_std = self.geneNormal.std(axis=1)
        geneCancer_std = self.geneCancer.std(axis=1)

        geneNormal_zfpkm = (geneNormal_log2.sub(geneNormal_mean, axis=0)).divide(geneNormal_std, axis=0)
        geneCancer_zfpkm = (geneCancer_log2.sub(geneCancer_mean, axis=0)).divide(geneCancer_std, axis=0)

        geneNormal_mean_zfpkm = geneNormal_zfpkm.mean(axis=1)
        geneCancer_mean_zfpkm = geneCancer_zfpkm.mean(axis=1)

        geneNormalIndex_less = geneNormal_mean_zfpkm[geneNormal_mean_zfpkm.values < (-3)].index
        geneCancerIndex_less = geneCancer_mean_zfpkm[geneCancer_mean_zfpkm.values < (-3)].index
        gene_less = geneNormalIndex_less.intersection(geneCancerIndex_less)

        self.geneNormal = self.geneNormal.drop(gene_less, axis=0, errors='ignore')
        self.geneCancer = self.geneCancer.drop(gene_less, axis=0, errors='ignore')

    def common_patients(self):
        Cancerpatients = pd.merge(pd.DataFrame(self.miRNACancer.columns), pd.DataFrame(self.geneCancer.columns))[0].tolist()
        Normalpatients = pd.merge(pd.DataFrame(self.miRNANormal.columns), pd.DataFrame(self.geneNormal.columns))[0].tolist()

        self.miRNACancer = self.miRNACancer[Cancerpatients]
        self.miRNANormal = self.miRNANormal[Normalpatients]
        self.geneCancer = self.geneCancer[Cancerpatients]
        self.geneNormal = self.geneNormal[Normalpatients]

    def matrices(self):
        Normal = pd.concat([self.geneNormal, self.miRNANormal], axis=0).sort_index().fillna(0)
        Cancer = pd.concat([self.geneCancer, self.miRNACancer], axis=0).sort_index().fillna(0)

        os.makedirs('./Matrices', exist_ok=True)
        Normal.to_csv(f'./Matrices/normal_{self.tumor}.csv', sep=';')
        Cancer.to_csv(f'./Matrices/cancer_{self.tumor}.csv', sep=';')

        Indices = pd.DataFrame(np.arange(0, len(Normal)), Normal.index)
        Indices.to_csv(f'./Matrices/indices_{self.tumor}.csv', sep=';')

        return Normal, Cancer, Indices
