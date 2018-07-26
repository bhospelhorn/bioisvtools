package waffleoRai_Utils;

import java.util.*;

/**
 * A template class for a simple binary tree structure.
 * @author Blythe Hospelhorn
 *
 * @param <T> Data type tree holds.
 * @version 1.0.0
 * @since March 2017
 */
public class BinTree<T> 
{
	
	class BinNode<Ty>
	{
		private Ty data;
		
		private BinNode<Ty> parent;
		private BinNode<Ty> right;
		private BinNode<Ty> left;
		
		BinNode(Ty newData)
		{
			this.data = newData;
			this.parent = null;
			this.right = null;
			this.left = null;
		}
		BinNode(Ty newData, BinNode<Ty> parent)
		{
			this.data = newData;
			this.parent = parent;
			this.left = null;
			this.right = null;
		}

		
		Ty getData()
		{
			return this.data;
		}
		void setData(Ty newData)
		{
			this.data = newData;
		}
		BinNode<Ty> getParent()
		{
			return this.parent;
		}
		public void setParent(BinNode<Ty> newParent)
		{
			this.parent = newParent;
		}
		public BinNode<Ty> getLeft()
		{
			return this.left;
		}
		public void setLeft(BinNode<Ty> newChild)
		{
			this.left = newChild;
		}
		public BinNode<Ty> getRight()
		{
			return this.right;
		}
		public void setRight(BinNode<Ty> newChild)
		{
			this.right = newChild;
		}

		public int numberChildren()
		{
			int count = 0;
			if (this.left != null) count++;
			if (this.right != null) count++;
			
			return count;
		}
		public boolean hasRightChild()
		{
			if (this.right != null) return true;
			return false;
		}
		public boolean hasLeftChild()
		{
			if (this.left != null) return true;
			return false;
		}
		public boolean hasParent()
		{
			if (this.parent != null) return true;
			return false;
		}
		
	}
	public static class BinNodeInfo<Ty>
	{
		private Ty data;
		private String HuffAddress;
		
		public BinNodeInfo(Ty theData, String theAddress)
		{
			this.data = theData;
			this.HuffAddress = theAddress;
		}
		
		public Ty getData()
		{
			return this.data;
		}
		public String getHuffCode()
		{
			return this.HuffAddress;
		}
		public void setData(Ty newData)
		{
			this.data = newData;
		}
		public void setHuff(String newAddress)
		{
			this.HuffAddress = newAddress;
		}
	}
	
	private BinNode<T> root;
	private BinNode<T> current;
	
	public BinTree()
	{
		/*Create empty tree*/
		root = null;
		current = null;
	}
	public BinTree(T rootData)
	{
		/*Create tree with nothing but a root.*/
		this.root = new BinNode<T>(rootData);
		this.current = this.root;
	}

	/*Insertion*/
	public void insertChildNode(T newData, boolean direction) /*Left = false, right = true*/
	{
		BinNode<T> child = new BinNode<T>(newData, this.current);
		
		if (direction)
		{
			if (!current.hasRightChild())
			{
				current.setRight(child);
				child.setParent(this.current);
			}
			else
			{
				child.setRight(current.getRight());
				child.getRight().setParent(child);
				current.setRight(child);
				child.setParent(this.current);
			}
		}
		else
		{
			if (!current.hasLeftChild())
			{
				current.setLeft(child);
			}
			else
			{
				child.setLeft(current.getLeft());
				child.getLeft().setParent(child);
				current.setLeft(child);
			}
		}
	}
	public void insertChildNode(BinTree<T> newBranch, boolean direction)
	{
		BinNode<T> branch = newBranch.getRootNode();
		BinNode<T> child = null;
		
		if (direction)
		{
			if (!current.hasRightChild())
			{
				current.setRight(branch);
				branch.setParent(this.current);
			}
			else
			{
				child = this.current.getRight();
				branch.setParent(this.current);
				this.current.setRight(branch);
				this.reattachBranch(child, true);
			}
		}
		else
		{
			if (!current.hasLeftChild())
			{
				current.setLeft(branch);
				branch.setParent(this.current);
			}
			else
			{
				child = this.current.getLeft();
				branch.setParent(this.current);
				this.current.setLeft(branch);
				this.reattachBranch(child, false);
			}
		}
	}
	public void insertParentNode(T newData, boolean direction)
	{
		BinNode<T> newParent = new BinNode<T>(newData);
		
		newParent.setParent(this.current.getParent());
		if (direction)
		{
			newParent.setRight(this.current);
			this.current.setParent(newParent);
		}
		else
		{
			newParent.setLeft(this.current);
			this.current.setParent(newParent);
		}
	}
	public void insertParentNode(BinTree<T> newBranch, boolean direction)
	{
		BinNode<T> branch = newBranch.getRootNode();
		BinNode<T> child = this.current;
		BinNode<T> oldParent = this.current.getParent();
		
		branch.setParent(oldParent);
		if (direction)
		{
			if (oldParent != null)
			{
				if (oldParent.getRight() == child) oldParent.setRight(branch);
				else oldParent.setLeft(branch);	
			}
			this.reattachBranch(child, true);
		}
		else
		{
			if (oldParent != null)
			{
				if (oldParent.getRight() == child) oldParent.setRight(branch);
				else oldParent.setLeft(branch);	
			}
			this.reattachBranch(child, false);
		}		
	}
	private void reattachBranch(BinNode<T> branchBase, boolean direction)
	{
		boolean found = false;
		BinNode<T> now = null;
		if (branchBase.hasParent())
		{
			/*First, start with parent and look for a place to reattach.*
			 */
			now = branchBase.getParent();	
		}
		else now = this.root;
		
		while (!found)
		{
			switch(now.numberChildren())
			{
				case 0:
					if (direction) now.setRight(branchBase);
					else now.setLeft(branchBase);
					branchBase.setParent(now);
					found = true;
					return;
				case 1:
					if (now.getRight() == null) now.setRight(branchBase);
					else now.setLeft(branchBase);
					branchBase.setParent(now);
					found = true;
					return;
				case 2:
					if (direction) now = now.getRight();
					else now = now.getLeft();
					break;
				default:
					break;
			}
		}
	}
	
	/*Deletion*/
	public void deleteCurrentNode(boolean heir)
	{
		BinNode<T> target = this.current;
		BinNode<T> p = target.getParent();
		BinNode<T> child = null;
		BinNode<T> sister = null;
		
		switch (target.numberChildren())
		{
		case 0:
			if (p.getRight() == target) p.setRight(null);
			else p.setLeft(null);
			this.current = p;
			return;
		case 1:
			if (target.hasLeftChild()) child = target.getLeft();
			else child = target.getRight();
			if (p.getRight() == target)
			{
				p.setRight(child);
				child.setParent(p);
			}
			else
			{
				p.setLeft(child);
				child.setParent(p);
			}
			this.current = child;
			return;
		case 2:
			if (heir)
			{
				child = target.getRight();
				sister = target.getLeft();
			}
			else
			{
				child = target.getLeft();
				sister = target.getRight();
			}
			if (p.getLeft() == target) p.setLeft(child);
			else p.setRight(child);
			child.setParent(p);
			sister.setParent(p);
			this.reattachBranch(sister, heir);
			this.current = child;
			return;
		default:
			return;
		}
	}
	public void deleteCurrentBranch()
	{
		BinNode<T> target = this.current;
		BinNode<T> p = target.getParent();
		
		if (!target.hasParent() && target == this.root)
		{
			this.root = null;
			return;
		}
		
		if (p.getLeft() == target) p.setLeft(null);
		else p.setRight(null);
		this.current = p;	
	}
	public void makeCurrentRoot()
	{
		this.root = this.current;
		this.current.setParent(null);
	}
	
	/*Movement*/
	public boolean hasRoot()
	{
		if (this.root != null) return true;
		return false;
	}
	public boolean moveLeft()
	{
		if (this.current.hasLeftChild())
		{
			this.current = this.current.getLeft();
			return true;
		}
		return false;
	}
	public boolean moveRight()
	{
		if (this.current.hasRightChild())
		{
			this.current = this.current.getRight();
			return true;
		}
		return false;	
	}
	public boolean moveUp()
	{
		if (this.current.hasParent())
		{
			this.current = this.current.getParent();
			return true;
		}
		return false;		
	}
	public boolean atRoot()
	{
		if (current == root) return true;
		return false;
	}
	public void moveToRoot()
	{
		this.current = this.root;
	}
	
	/*Information retrieval and management*/
	public T getCurrentData()
	{
		return this.current.getData();
	}
	public boolean currentIsLeaf()
	{
		if (current.numberChildren() > 0) return false;
		return true;
	}
	public int currentChildren()
	{
		return this.current.numberChildren();
	}
	public void setCurrentData(T newData)
	{
		this.current.setData(newData);
	}
	public T getRootData()
	{
		return this.root.getData();
	}
	public int rootChildren()
	{
		return this.root.numberChildren();
	}
	public void setRootData(T newData)
	{
		this.root.setData(newData);
	}

	private boolean searcher(T target, BinNode<T> parent)
	{		
		/*Check self*/
		if (parent.getData().equals(target))
		{
			this.current = parent;
			return true;
		}
			
		/*Check left branch*/
		if (parent.hasLeftChild())
		{
			if (searcher(target, parent.getLeft()))
			{
				return true;
			}
		}
		
		/*Check right branch*/
		if (parent.hasRightChild())
		{
			if (searcher(target, parent.getRight()))
			{
				return true;
			}
		}
		
		return false;
	}
	
	public boolean searchAndSetCurrent(T target)
	{
		return this.searcher(target, this.root);
	}
	public String searchAndGetLocation(T target)
	{
		String hadd = "NOT FOUND";
		BinNode<T> temp = this.current;
		if (searchAndSetCurrent(target))
		{
			hadd = huffAddress(this.current);
		}
		
		this.current = temp;
		
		return hadd;
	}
	
	BinNode<T> getRootNode()
	{
		return this.root;
	}
	BinNode<T> getCurrentNode()
	{
		return this.current;
	}
	
	private void addNodeToList(List<BinNodeInfo<T>> myList, BinNode<T> myNode)
	{
		if (myNode.hasLeftChild()) addNodeToList(myList, myNode.getLeft());
		BinNodeInfo<T> myInfo = new BinNodeInfo<T>(myNode.getData(), huffAddress(myNode));
		myList.add(myInfo);
		if (myNode.hasRightChild()) addNodeToList(myList, myNode.getRight());
	}
	public List<BinNodeInfo<T>> toList()
	{
		/*Returns a representation of the tree.*/
		List<BinNodeInfo<T>> myList = new ArrayList<BinNodeInfo<T>>();
		addNodeToList(myList, this.root);
		
		return myList;
	}
	
	private String huffAddress(BinNode<T> node)
	{
		/*Returns a Huffman style binary representation of the node's location*/
		String hLocation = "";
		BinNode<T> now = node;
		BinNode<T> last = null;
		
		while (now != root)
		{
			last = now;
			now = now.getParent();
			if (last == now.getLeft()) hLocation = "0" + hLocation;
			if (last == now.getRight()) hLocation = "1" + hLocation;
		}
		
		return hLocation;
	}
	
	public String currentString()
	{
		/*Returns a string representing the current node.*/
		if (this.current == null) return "Current node is null!\n";
		
		String sdata = "Data: " + this.current.getData().toString() + "\n";
		String sparent = null;
		String sleft = null;
		String sright = null;
		if (this.current.hasParent())
		{
			sparent = "Parent: (" + this.current.getParent().getData().toString() + ")\n";
		}
		else sparent = "Parent: None\n";
		if (this.current.hasLeftChild())
		{
			sleft = "Left Child: (" + this.current.getLeft().getData().toString() + ")\n";
		}
		else sleft = "Left Child: None\n";
		if (this.current.hasRightChild())
		{
			sright = "Right Child: (" + this.current.getRight().getData().toString() + ")\n";
		}
		else sright = "Right Child: None\n";
			
		
		return sdata + sparent + sleft + sright;
	}
	public String currentLocation()
	{
		/*Returns a Huffman style binary representation of the current node's location*/
		String hLocation = huffAddress(this.current);
		
		return hLocation;
	}
	private String branchString(BinNode<T> branchBase)
	{
		String s = "";
		
		if (branchBase.hasLeftChild()) s = s + branchString(branchBase.getLeft());
		s = s + "TREE ADDRESS: " + huffAddress(branchBase) + "\n" + branchBase.getData().toString() + "\n";
		if (branchBase.hasRightChild()) s = s + branchString(branchBase.getRight());
		
		return s;
	}
	public String toString()
	{
		return this.branchString(this.root);
	}
	
}
