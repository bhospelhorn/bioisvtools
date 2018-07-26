package hospelhornbg_genomeBuild;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.UCSCGVBED;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_svtools.ConsoleMain;

public class MillsConvert {
	
	private static class Record
	{
		public String sample;
		public String chrom;
		public int start;
		public int end;
		public String type;
		public String source;
		
		public String toString()
		{
			return (chrom + ":" + start + "-" + end + " [" + type + "]");
		}
	}

	public static void main(String[] args) {
		
		//String[] iargs = {"-d", "C:\\Users\\hospelhornbg\\Documents\\bioisvtools"};
		String homedir = ConsoleMain.getUserHome(true);
		//Install.installBIOISVTOOLS(iargs, homedir);
		
		//Load the genome builds...
		GenomeBuild hg18 = ConsoleMain.loadBuild(homedir, "hg18", true);
		GenomeBuild hg19 = ConsoleMain.loadBuild(homedir, "hg19", true);
		GenomeBuild hg38 = ConsoleMain.loadBuild(homedir, "hg38", true);
		
		String delpath = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\nature09708-s6_DEL.csv";
		String duppath = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\nature09708-s6_DUPINS.csv";
		
		LinkedList<Record> records = new LinkedList<Record>();
		
		//Read in raw files!
		try
		{
			FileReader fr = new FileReader(delpath);
			BufferedReader br = new BufferedReader(fr);
			int linecount = 0;
			String line = br.readLine(); //Skip header line
			while ((line = br.readLine()) != null)
			{
				linecount++;
				String[] fields = line.split(",");
				if (fields == null || fields.length < 1) {
					System.err.println("ERROR (del): Line " + linecount + " could not be read!");
					continue;
				}
				if (fields.length != 6) {
					System.err.println("ERROR (del): Line " + linecount + " could not be parsed!");
					continue;
				}
				Record r = new Record();
				r.sample = fields[0];
				r.chrom = fields[1];
				r.type = fields[4];
				r.source = fields[5];
				try
				{
					r.start = Integer.parseInt(fields[2]);
					r.end = Integer.parseInt(fields[3]);
				}
				catch (NumberFormatException e)
				{
					System.err.println("ERROR (del): Line " + linecount + " could not be parsed (Number error)!");
					continue;
				}
				records.add(r);
				
			}
			br.close();
			fr.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
			System.exit(1);
		}
		
		//Repeat with dup.
		//Copypaste because I'm too lazy to wire up a proper function
		try
		{
			FileReader fr = new FileReader(duppath);
			BufferedReader br = new BufferedReader(fr);
			int linecount = 0;
			String line = br.readLine(); //Skip header line
			while ((line = br.readLine()) != null)
			{
				linecount++;
				String[] fields = line.split(",");
				if (fields == null || fields.length < 1) {
					System.err.println("ERROR (dup/ins): Line " + linecount + " could not be read!");
					continue;
				}
				if (fields.length != 6) {
					System.err.println("ERROR (dup/ins): Line " + linecount + " could not be parsed!");
					continue;
				}
				Record r = new Record();
				r.sample = fields[0];
				r.chrom = fields[1];
				r.type = fields[4];
				r.source = fields[5];
				try
				{
					r.start = Integer.parseInt(fields[2]);
					r.end = Integer.parseInt(fields[3]);
				}
				catch (NumberFormatException e)
				{
					System.err.println("ERROR (dup/ins): Line " + linecount + " could not be parsed (Number error)!");
					continue;
				}
				records.add(r);
				
			}
			br.close();
			fr.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
			System.exit(1);
		}
		
		//Prepare variant pools
		VariantPool vp18 = new VariantPool(1);
		vp18.setGenomeBuild(hg18);
		VariantPool vp19 = new VariantPool(1);
		vp19.setGenomeBuild(hg19);
		VariantPool vp38 = new VariantPool(1);
		vp38.setGenomeBuild(hg38);
		VariantPool[] pools = {vp18, vp19, vp38};
		
		InfoDefinition sourcedef = new InfoDefinition("1000GSOURCE", VariantPool.INFODEF_STRING, "The source of this variant call as listed in the spreadsheet provided by Mills et. al.", 1);
		for (VariantPool pool : pools) {
			pool.addSample("NA12878"); //GIAB
			//NA12156 is included too, but she appears to be some rando?
			StructuralVariant.addStandardDefs(pool, false);
			pool.addInfoField(sourcedef.getKey(), sourcedef);
			pool.addFormatField(Genotype.INFODEF_GT.getKey(), Genotype.INFODEF_GT);
		}
		
		//Now we try to turn these records into variants...
		int count = 0;
		for (Record r : records)
		{
			count++;
			if (!r.sample.equals("NA12878")) continue;
			if (r.chrom.equals("chr23")) r.chrom = "chrX";
			else if (r.chrom.equals("23")) r.chrom = "X";
			SVType t = SVType.getType(r.type);
			if (t == null)
			{
				System.err.println("Type (" + r.type + ") for variant " + r.chrom + ":" + r.start + "-" + r.end + " invalid! Passing over...");
				continue;
			}
			//The fun part... Checking the contigs!
			Contig c18 = hg18.getContig(r.chrom);
			Contig c19 = hg19.getContig(r.chrom);
			Contig c38 = hg38.getContig(r.chrom);
			if (c18 == null)
			{
				System.err.println("ERROR: contig " + r.chrom + " does not appear to exist in hg18!");
				System.err.println("Tossing variant " + r.toString() + " from hg18...");
			}
			if (c19 == null)
			{
				System.err.println("ERROR: contig " + r.chrom + " does not appear to exist in hg19!");
				System.err.println("Tossing variant " + r.toString() + " from hg19...");
			}
			if (c38 == null)
			{
				System.err.println("ERROR: contig " + r.chrom + " does not appear to exist in hg38!");
				System.err.println("Tossing variant " + r.toString() + " from hg38...");
			}
			
			if (c18 != null)
			{
				boolean s18 = (r.start < c18.getLength());
				boolean e18 = (r.end <= c18.getLength());
				if (!s18)
				{
					System.err.println("Variant rejected for build hg18: " + r.toString());
					System.err.println("Reason: Start " + r.start + " exceeds length for contig " + c18.getUDPName() + " (" + c18.getLength() + ")");
				}
				if (!e18)
				{
					System.err.println("Variant rejected for build hg18: " + r.toString());
					System.err.println("Reason: End " + r.end + " exceeds length for contig " + c18.getUDPName() + " (" + c18.getLength() + ")");
				}
				if (s18 && e18)
				{
					StructuralVariant sv = new StructuralVariant();
					sv.setType(t);
					sv.setChromosome(c18);
					sv.setPosition(r.start);
					sv.setEndPosition(r.end);
					sv.setRefAllele("N");
					sv.addAltAllele("<" + t.toString() + ">");
					sv.addGenotype("NA12878", new Genotype());
					int len = r.end - r.start;
					if (t == SVType.DEL) len *= -1;
					sv.setSVLength(len);
					sv.addInfoField(r.source, sourcedef);
					sv.setVariantName(t.toString() + "_" + count);
					vp18.addVariant(sv);
				}
			}
			if (c19 != null)
			{
				boolean s19 = (r.start < c19.getLength());
				boolean e19 = (r.end <= c19.getLength());
				if (!s19)
				{
					System.err.println("Variant rejected for build hg19: " + r.toString());
					System.err.println("Reason: Start " + r.start + " exceeds length for contig " + c19.getUDPName() + " (" + c19.getLength() + ")");
				}
				if (!e19)
				{
					System.err.println("Variant rejected for build hg19: " + r.toString());
					System.err.println("Reason: End " + r.end + " exceeds length for contig " + c19.getUDPName() + " (" + c19.getLength() + ")");
				}
				if (s19 && e19)
				{
					StructuralVariant sv = new StructuralVariant();
					sv.setType(t);
					sv.setChromosome(c19);
					sv.setPosition(r.start);
					sv.setEndPosition(r.end);
					sv.setRefAllele("N");
					sv.addAltAllele("<" + t.toString() + ">");
					sv.addGenotype("NA12878", new Genotype());
					int len = r.end - r.start;
					if (t == SVType.DEL) len *= -1;
					sv.setSVLength(len);
					sv.addInfoField(r.source, sourcedef);
					sv.setVariantName(t.toString() + "_" + count);
					vp19.addVariant(sv);
				}
			}
			if (c38 != null)
			{
				boolean s38 = (r.start < c38.getLength());
				boolean e38 = (r.end <= c38.getLength());
				if (!s38)
				{
					System.err.println("Variant rejected for build hg38: " + r.toString());
					System.err.println("Reason: Start " + r.start + " exceeds length for contig " + c38.getUDPName() + " (" + c38.getLength() + ")");
				}
				if (!e38)
				{
					System.err.println("Variant rejected for build hg38: " + r.toString());
					System.err.println("Reason: End " + r.end + " exceeds length for contig " + c38.getUDPName() + " (" + c38.getLength() + ")");
				}
				if (c38 != null && s38 && e38)
				{
					StructuralVariant sv = new StructuralVariant();
					sv.setType(t);
					sv.setChromosome(c38);
					sv.setPosition(r.start);
					sv.setEndPosition(r.end);
					sv.setRefAllele("N");
					sv.addAltAllele("<" + t.toString() + ">");
					sv.addGenotype("NA12878", new Genotype());
					int len = r.end - r.start;
					if (t == SVType.DEL) len *= -1;
					sv.setSVLength(len);
					sv.addInfoField(r.source, sourcedef);
					sv.setVariantName(t.toString() + "_" + count);
					vp38.addVariant(sv);
				}
			}
			
			
		}
		
		//Output all three pools as vcfs and BED tracks
		String vcf18 = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\MillsConvert\\NA12878_1000G_as_hg18.vcf";
		String vcf19 = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\MillsConvert\\NA12878_1000G_as_hg19.vcf";
		String vcf38 = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\MillsConvert\\NA12878_1000G_as_hg38.vcf";
		
		String bed18 = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\MillsConvert\\NA12878_1000G_as_hg18.bed";
		String bed19 = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\MillsConvert\\NA12878_1000G_as_hg19.bed";
		String bed38 = "X:\\usr\\hospelhornbg\\GIAB\\NA12878\\MillsConvert\\NA12878_1000G_as_hg38.bed";
		
		try 
		{
			vp18.sortVariants();
			VCF.writeVCF(vp18, "[NHGRIUDP]", vcf18);
			vp19.sortVariants();
			VCF.writeVCF(vp19, "[NHGRIUDP]", vcf19);
			vp38.sortVariants();
			VCF.writeVCF(vp38, "[NHGRIUDP]", vcf38);
			
			UCSCGVBED bed = new UCSCGVBED("NA12878 Mills et al. i_hg18", vp18.getVariants());
			bed.setDescription("Mills et al. NA12878 SV callset interpreted using build hg18");
			bed.write(bed18, "NA12878");
			
			bed = new UCSCGVBED("NA12878 Mills et al. i_hg19", vp19.getVariants());
			bed.setDescription("Mills et al. NA12878 SV callset interpreted using build hg19");
			bed.write(bed19, "NA12878");
			
			bed = new UCSCGVBED("NA12878 Mills et al. i_hg38", vp38.getVariants());
			bed.setDescription("Mills et al. NA12878 SV callset interpreted using build hg38");
			bed.write(bed38, "NA12878");
			
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		
		
		//Print contig info...
		List<Contig> ctg18 = hg18.getChromosomes();
		System.out.println("hg18 contigs...");
		for (Contig c : ctg18)
		{
			System.out.println(c.printInfo());
		}

	}

}
