package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.FileReader;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;

public class GIABTruthTest {

	public static void main(String[] args) {
		
		String infile = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\MillsConvert\\NA12878_1000G_hg18_to_hg19_liftover.bed";
		String outfile = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\MillsConvert\\NA12878_1000G_hg18_to_hg19_liftover.vcf";
		
		VariantPool pool = new VariantPool(1);
		pool.addSample("NA12878");
		InfoDefinition def = StructuralVariant.INFODEF_INFO_SVTYPE;
		pool.addInfoField(def.getKey(), def);
		def = StructuralVariant.INFODEF_INFO_SVLEN;
		pool.addInfoField(def.getKey(), def);
		def = StructuralVariant.INFODEF_INFO_END;
		pool.addInfoField(def.getKey(), def);
		def = Genotype.INFODEF_GT;
		pool.addFormatField(def.getKey(), def);
		String alt = StructuralVariant.INFODEF_ALT_DEL;
		pool.addCustomAlt(alt, StructuralVariant.getAltDefDescription(alt));
		alt = StructuralVariant.INFODEF_ALT_DUP;
		pool.addCustomAlt(alt, StructuralVariant.getAltDefDescription(alt));
		alt = StructuralVariant.INFODEF_ALT_INS;
		pool.addCustomAlt(alt, StructuralVariant.getAltDefDescription(alt));
		
		//Load Genome Build
		String gbPath = "C:\\Users\\hospelhornbg\\Desktop\\GRCh37.gbdh";
		
		try 
		{
			GenomeBuild gb = new GenomeBuild(gbPath);
			pool.setGenomeBuild(gb);
			
			FileReader fr = new FileReader(infile);
			BufferedReader br = new BufferedReader(fr);
			
			//Skip the first line
			br.readLine();
			String line = null;
			while ((line = br.readLine()) != null)
			{
				StructuralVariant sv = new StructuralVariant();
				sv.setRefAllele("N");
				//Tab split 
				String[] fields = line.split("\t");
				Contig c = gb.getContig(fields[0]);
				int pos = Integer.parseInt(fields[1]);
				int end = Integer.parseInt(fields[2]);
				String varname = fields[3];
				//Figure out type
				SVType t = SVType.DEL;
				if (varname.contains("DUP")) t = SVType.DUP;
				if (varname.contains("INS")) t = SVType.INS;
				//Calculate length
				int len = end - pos;
				if (t == SVType.DEL) len = pos - end;
				sv.setChromosome(c);
				sv.setEndPosition(end);
				sv.setPosition(pos);
				sv.setVariantName(varname);
				sv.addAltAllele("<" + t.toString() + ">");
				sv.setType(t);
				sv.setSVLength(len);
				sv.addGenotype("NA12878", new Genotype());
				pool.addVariant(sv);
			}
			
			br.close();
			fr.close();
			
			pool.sortVariants();
			VCF.writeVCF(pool, "Bioisvtools[NIHUDP]", outfile);
			
			
		} 
		catch (Exception e) 
		{
			e.printStackTrace();
		}
	}

}
