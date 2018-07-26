package hospelhornbg_bioinformatics;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Random;

public class UCSCGVBED {
	
	//chrom	start	end	name	score	strand	thickstart	thickend	RGB	blockCount	blockSizes	blockStarts
	
	public static final String COLOR_DEL = "231,16,26"; //Red
	public static final String COLOR_DELME = "129,10,16"; //Dark red
	public static final String COLOR_DUP = "43,190,28"; //Green
	public static final String COLOR_DUPTANDEM = "157,217,27"; //Spring green
	public static final String COLOR_INS = "29,65,198"; //Blue
	public static final String COLOR_INSME = "11,33,115"; //Dark blue
	public static final String COLOR_INV = "180,63,144"; //Purple
	public static final String COLOR_BND = "248,248,20"; //Yellow
	public static final String COLOR_TRA = "234,111,10"; //Orange
	public static final String COLOR_SNV = "0,0,0"; //Black
	public static final String COLOR_CNVUNK = "92,24,95"; //Dark purple
	
	private Collection<Variant> variants;
	private Collection<String> sampleList; //For bookkeeping, mostly. Optional.
	
	private String name;
	private String description;
	
	
	public UCSCGVBED(String trackname, Collection<Variant> pool)
	{
		variants = new ArrayList<Variant>(pool.size());
		variants.addAll(pool);
		name = trackname;
		if (name == null || name.isEmpty()) genName();
		description = "";
	}
	
	private void genName()
	{
		Random rand = new Random();
		int i = rand.nextInt();
		name = String.format("Track%08x", i);
	}
	
	public String getName()
	{
		return name;
	}
	
	public String getDescription()
	{
		return description;
	}
	
	public Collection<String> getSampleList()
	{
		return sampleList;
	}
	
	public void setName(String newName)
	{
		if (newName == null) return;
		if (newName.isEmpty()) return;
		name = newName;
	}
	
	public void setDescription(String desc)
	{
		if (desc == null) return;
		description = desc;
	}
	
	public void setVariants(Collection<Variant> pool)
	{
		if (pool == null || pool.isEmpty()) variants.clear();
		variants = new ArrayList<Variant>(pool.size());
		variants.addAll(pool);
	}
	
	public void setSampleList(Collection<String> samples)
	{
		sampleList = samples;
	}
	
	public static String getSVColor(SVType t)
	{
		switch(t)
		{
		case BED_REGION:
			return COLOR_SNV;
		case BND:
			return COLOR_BND;
		case CNV:
			return COLOR_CNVUNK;
		case DEL:
			return COLOR_DEL;
		case DELME:
			return COLOR_DELME;
		case DUP:
			return COLOR_DUP;
		case INS:
			return COLOR_INS;
		case INSME:
			return COLOR_INSME;
		case INV:
			return COLOR_INV;
		case OTHER:
			return COLOR_SNV;
		case TANDEM:
			return COLOR_DUPTANDEM;
		default:
			return COLOR_SNV;
		}
	}
	
	private String writeHeader()
	{
		//track name=[name] description=[desc] visibility=n itemRgb="On"
		return "track name=\"" + this.name + "\" description=\"" + this.description + "\" visibility=2 itemRgb=\"On\"";
	}
	
	public void write(String filepath, String sampleName) throws IOException
	{
		if (filepath == null || filepath.isEmpty()) throw new IOException();
		if (variants == null) return;
		FileWriter writer = new FileWriter(filepath);
		writer.write(this.writeHeader() + "\n");
		int counter = 0;
		for (Variant v : variants)
		{
			String line = v.toViewerBEDLine(sampleName);
			if (line != null) {
				if (counter < variants.size() - 1) writer.write(line + "\n");
				else writer.write(line);	
			}
			counter++;
		}
		
		
		writer.close();
	}
	
	
}
