"""
Excel Data Reader Module

Handles reading and processing Excel files for the ACMG SF Variants Report.
Maps Excel data to template variables according to the data mapping guide.
"""

import pandas as pd
from datetime import datetime
from typing import Dict, List, Any, Optional


class ExcelDataReader:
    """Reads and processes Excel data for report generation"""

    def __init__(self, excel_path: str, debug: bool = False):
        self.excel_path = excel_path
        self.debug = debug
        self.required_tabs = {
            'Sample Summary': True,
            'ACMG SF (P-LP)': True,
            'Coverage gaps': False,  # Optional
            'ACMG Genes Coverage': False,  # Optional
        }

    def read_all_tabs(self) -> Dict[str, pd.DataFrame]:
        """Read all relevant tabs from Excel file"""
        if self.debug:
            print(f"Reading Excel file: {self.excel_path}")

        try:
            # Read the Excel file
            excel_data = pd.read_excel(self.excel_path, sheet_name=None)

            if self.debug:
                print(f"Available tabs: {list(excel_data.keys())}")

            # Check for required tabs
            for tab_name, required in self.required_tabs.items():
                if tab_name not in excel_data and required:
                    raise ValueError(f"Required tab '{tab_name}' not found in Excel file")

            return excel_data

        except Exception as e:
            raise Exception(f"Error reading Excel file: {e}")

    def validate_data(self, data: Dict[str, pd.DataFrame]) -> None:
        """Validate that required columns exist in each tab"""
        validations = {
            'Sample Summary': ['Sample_ID'],  # Add other required columns
            'ACMG SF (P-LP)': ['Gene', 'HGVSc', 'HGVSp', 'ClinVar'],
        }

        for tab_name, required_columns in validations.items():
            if tab_name in data:
                df = data[tab_name]
                missing_columns = [col for col in required_columns if col not in df.columns]
                if missing_columns:
                    raise ValueError(f"Missing required columns in '{tab_name}': {missing_columns}")

    def process_for_template(self, data: Dict[str, pd.DataFrame], sample_id: str) -> Dict[str, Any]:
        """Process Excel data into format needed for HTML template"""

        template_data = {
            'sample_id': sample_id,
            'report_generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
        }

        # Process Page 1 header data
        template_data.update(self._process_header_data(data))

        # Process variant data for Page 1 & 2
        template_data.update(self._process_variant_data(data))

        # Process QC metrics for Page 3
        template_data.update(self._process_qc_metrics(data))

        # Process coverage data for Page 3
        template_data.update(self._process_coverage_data(data))

        return template_data

    def _process_header_data(self, data: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
        """Extract header information from Sample Summary tab"""
        if 'Sample Summary' not in data:
            return {}

        df = data['Sample Summary']
        header_data = {}

        # Map Excel columns to template variables based on actual lean report structure
        column_mapping = {
            'Sample_ID': 'sample_id',
            'Assay': 'assay',
            'Build': 'reference_build',
        }

        # Extract first row of data (assuming single sample per file)
        if len(df) > 0:
            for excel_col, template_var in column_mapping.items():
                if excel_col in df.columns:
                    header_data[template_var] = df[excel_col].iloc[0]

        # Set defaults for missing data
        header_data.setdefault('assay', 'WES')
        header_data.setdefault('reference_build', 'GRCh38')
        header_data.setdefault('pipeline', 'Bionl_Lean_call v1.0')
        header_data.setdefault('databases', 'ClinVar (2025-01), gnomAD v4.1, Ensembl VEP Release 115, REVEL (latest release), AlphaMissense (Science 2023, updated 2025-05)')

        return header_data

    def _process_variant_data(self, data: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
        """Extract and filter variant data for Page 1 and 2"""
        if 'ACMG SF (P-LP)' not in data:
            return {
                'page1_variants': [],
                'page2_variants': []
            }

        df = data['ACMG SF (P-LP)']
        # Filter for pathogenic/likely pathogenic variants
        pathogenic_mask = df['ClinVar'].isin(['Pathogenic', 'Likely pathogenic','Conflicting_classifications_of_pathogenicity'])
        pathogenic_variants = df[pathogenic_mask]

        # Split based on star rating (>=2 stars for Page 1, <2 stars for Page 2)
        page1_variants = []
        page2_variants = []

        for _, row in pathogenic_variants.iterrows():
            variant_data = {
                'gene': row.get('Gene', '—'),
                'hgvsc': row.get('HGVSc', '—'),
                'hgvsp': row.get('HGVSp', '—'),
                'clinvar_id': row.get('ClinVar_Link', ''),  # This contains the ClinVar ID
            }

            # Check star rating
            stars = row.get('ClinVar_Stars', 0)
            try:
                stars = int(stars) if pd.notna(stars) else 0
            except (ValueError, TypeError):
                stars = 0

            if stars >= 2:
                page1_variants.append(variant_data)
            else:
                page2_variants.append(variant_data)

        return {
            'page1_variants': page1_variants,
            'page2_variants': page2_variants
        }

    def _process_qc_metrics(self, data: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
        """Extract QC metrics from Sample Summary tab"""
        if 'Sample Summary' not in data:
            return {}

        df = data['Sample Summary']
        print(df)
        metrics = {}

        # Map Excel columns to template variables based on actual lean report structure
        if len(df) > 0:
            row = df.iloc[0]
            
            # Handle coverage metrics
            mean_cov = row.get('Mean_Coverage', '—')
            if mean_cov != '—' and pd.notna(mean_cov):
                metrics['mean_coverage'] = f"{mean_cov:.2f}x"
            else:
                metrics['mean_coverage'] = '—'
            
            # Handle mapping percentage
            mapped_pct = row.get('Mapped_Percent', '—')
            if mapped_pct != '—' and pd.notna(mapped_pct):
                metrics['mapped_percent'] = f"{mapped_pct:.2f}%"
            else:
                metrics['mapped_percent'] = '—'
            
            # Handle duplicate percentage
            dup_pct = row.get('Duplicate_Percent', '—')
            if dup_pct != '—' and pd.notna(dup_pct):
                metrics['duplicate_percent'] = f"{dup_pct:.2f}%"
            else:
                metrics['duplicate_percent'] = '—'
            
            # Handle insert size
            insert_size = row.get('Mean_Insert_Size', '—')
            if insert_size != '—' and pd.notna(insert_size):
                metrics['mean_insert_size'] = f"{insert_size:.1f}bp"
            else:
                metrics['mean_insert_size'] = '—'
            
            # Handle ACMG coverage - this might be in a different sheet or calculated
            # For now, set a default or try to find it in coverage data
            metrics['acmg_covered_percent'] = row.get('ACMG_PctRegionsCovered', '—')  # Will be updated from coverage data if available

        return metrics

    def _process_coverage_data(self, data: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
        """Extract coverage data for both coverage gaps and overall coverage sections"""
        coverage_data = {
            'coverage_gaps': [],
            'overall_coverage': []
        }

        # Process coverage gaps from "Coverage gaps" sheet
        if 'Coverage gaps' in data:
            gaps_df = data['Coverage gaps']
            for _, row in gaps_df.iterrows():
                coverage_pct = row.get('Pct>=20x', 100)
                try:
                    coverage_pct = float(coverage_pct) if pd.notna(coverage_pct) else 100
                except (ValueError, TypeError):
                    coverage_pct = 100

                gene_data = {
                    'gene': row.get('Gene', '—'),
                    'transcript': row.get('MANE_ID', '—'),
                    'coverage': f"{coverage_pct:.1f}%"
                }

                # Add to gaps if coverage < 100%
                if coverage_pct < 100:
                    coverage_data['coverage_gaps'].append(gene_data)

        # Process overall coverage from "ACMG SF Genes Coverage" sheet
        if 'ACMG SF Genes Coverage' not in data:
            return coverage_data

        df = data['ACMG SF Genes Coverage']

        # Process overall coverage (length-weighted mean coverage per gene)
        for _, row in df.iterrows():
            coverage_pct = row.get('Pct20', 100)  # This is the length-weighted mean
            try:
                coverage_pct = float(coverage_pct) if pd.notna(coverage_pct) else 100
            except (ValueError, TypeError):
                coverage_pct = 100

            gene_data = {
                'gene': row.get('Gene', '—'),
                'transcript': row.get('MANE_ID', '—'),
                'coverage': f"{coverage_pct:.1f}%"
            }

            # Add all genes to overall coverage
            coverage_data['overall_coverage'].append(gene_data)

        return coverage_data