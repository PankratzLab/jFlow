package org.genvisis.one.george;

import static org.junit.Assert.*;

import org.junit.Test;

public class TreeTest {

	public static void printArray(Object[] arr) {
		for (int i = 0; i < arr.length; i++)
			if (i == arr.length-1)
				System.out.println(arr[i].toString());
			else
				System.out.print(arr[i].toString() + ",");
	}
	@Test
	public void testbase() {
		Node n = new Node(1, null);
		Tree t = new Tree(n, 1);
		try {
			Integer[] arr = t.getArray(1);
			printArray(arr);
		}
		catch (Exception e) {
			e.printStackTrace();
			System.out.println("exception");
			fail();
		}
		
	}
	
	@Test
	public void test2() {
		Node child1 = new Node(2, new Node[]{new Node(4, null)});
		Node child2 = new Node(3, new Node[]{new Node(6, null)});
		Node n = new Node(1, new Node[]{child1, child2});
		Tree t = new Tree(n, 6);
		try {
			Integer[] arr = t.getArray(6);
			printArray(arr);
		}
		catch (Exception e) {
			e.printStackTrace();
			System.out.println("exception");
			fail();
		}
	}

}
