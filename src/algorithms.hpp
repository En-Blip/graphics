
#ifndef ALGORITHMS_HPP_GRAPHICS

class TreeNode {
public:
    TreeNode(int value);
    ~TreeNode();
    int value;
    char balance;
    TreeNode* left;
    TreeNode* right;
    TreeNode* parent;
};

class AVLTree {
public:
    AVLTree();
    ~AVLTree();
    TreeNode* root;
    void insert(TreeNode* node);

private:
    void balanceSingleRotation(TreeNode* pivot, TreeNode* son, bool dir);
    void balanceDoubleRotation(TreeNode* pivot, TreeNode* son, TreeNode* grandson, bool dir);
    void deleteNodeRecursive(TreeNode* node);
};

class MaxHeap {
public:
    MaxHeap(int capacity);
    ~MaxHeap();
    void insert(int value);
    int* heapsort();
    int size;
    int capacity;
    int* heap;

private:
    void heapify(int from);
};

void print(TreeNode* root, int depth);

#define ALGORITHMS_HPP_GRAPHICS
#endif
