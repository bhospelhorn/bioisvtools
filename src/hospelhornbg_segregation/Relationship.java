package hospelhornbg_segregation;

public interface Relationship {
	
	public static final String[] DISTANCE_ENGLISH = {"First", "Second", "Third", "Fourth",
													 "Fifth", "Sixth", "Seventh", "Eighth",
													 "Ninth", "Tenth", "Eleventh", "Twelfth",
													 "Thirteenth", "Fourteenth", "Fifteenth", "Sixteenth",
													 "Seventeenth", "Eighteenth", "Nineteenth", "Twentieth"};
	
	public static final String[] REMOVAL_ENGLISH = {"Once", "Twice", "Thrice", "Four Times",
													"Five Times", "Six Times", "Seven Times", "Eight Times",
													"Nine Times", "Ten Times", "Eleven Times", "Twelve Times",
													"Thirteen Times", "Fourteen Times", "Fifteen Times", "Sixteen Times",
													"Seventeen Times", "Eighteen Times", "Nineteen Times", "Twenty Times"};
	
	public Individual getTarget();
	public Individual getRelative();
	public Individual getCommonAncestor();
	public Individual[] getCommonAncestors();
	
	public boolean isAncestor();
	public boolean isDescendant();
	public boolean isParent();
	public boolean isChild();
	public boolean isSibling();
	public boolean isHalfSibling();
	public boolean isSelf();
	
	public boolean commonAncestorCouple();
	
	public int targetGenerationsToCommonAncestor();
	public int relativeGenerationsToCommonAncestor();
	
	public String toString_English();

}
