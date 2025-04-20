#include "algorithms.hpp"
#include <iostream>
#include <cmath>

TreeNode::TreeNode(int value){
    this->value = value;
    this->parent = NULL;
    this->left = NULL;
    this->right = NULL;
}

TreeNode::~TreeNode(){
    if(this->parent != NULL){
        delete this->parent;
    }
    if(this->left != NULL){
        delete this->left;
    }
    if(this->right != NULL){
        delete this->right;
    }
}

AVLTree::AVLTree(){
   this->root = NULL; 
}

AVLTree::~AVLTree(){
    deleteNodeRecursive(this->root);
}

void AVLTree::deleteNodeRecursive(TreeNode* node){
    if(node == NULL) return;

    deleteNodeRecursive(node->left);
    deleteNodeRecursive(node->right);
    delete node;
}

void AVLTree::insert(TreeNode* node){
    node->balance = 0;
    if(root == NULL){
        root = node;
        return;
    }

    TreeNode* pivot = NULL;
    TreeNode* curNode = root;
    TreeNode* nextNode = NULL;
    do{
        if(nextNode != NULL) curNode = nextNode;
        if(node->value > curNode->value){
            nextNode = curNode->right;

        }else if(node->value <= curNode->value){
            nextNode = curNode->left;
            
        }

    }while(nextNode != NULL);

    nextNode = curNode;
    if(pivot == NULL){
        if(node->value > curNode ->value){
            curNode->right = node;
            node->parent = curNode;
        }else{
            curNode->left = node;
            node->parent = curNode;
        }

        while(nextNode != NULL){
            if(node->value > nextNode->value){
                nextNode->balance++;
            }else{
                nextNode->balance--;
            }
    
            if(pivot != NULL && abs(nextNode->balance) > 1){
                pivot = nextNode;
            }
    
            if(nextNode->balance == 0){
                break;
            }
            nextNode = nextNode->parent;
        }
    }else if((node->value > pivot->value && node->value < pivot->right->value) || (node->value < pivot->value && node->value > pivot->left->value)){
        bool dir = (node->value > pivot->value);

        if(dir){
            curNode->right = node;
            node->parent = curNode;
        }else{
            curNode->left = node;
            node->parent = curNode;
        }
        TreeNode* son = (dir)?pivot->right:pivot->left;
        TreeNode* grandson = (dir)?pivot->right->left:pivot->left->right;
        balanceDoubleRotation(pivot, son, grandson, dir);

        while(nextNode != son && nextNode != grandson){
           nextNode->balance += (node->value > nextNode->value)?1:-1; 
        }
        if(node->value > pivot->value){// rl rotation
            if(node->value > grandson->value){// node > grandsom
                pivot->balance = -1;
            }else{
                pivot->balance = 0;
                son->balance = 1;
            }
        }else{ // lr rotation
            if(node->value <= grandson->value){
                pivot->balance = -1;
            }else{
                pivot->balance = 0;
                son->balance = 1;
            }
        }
    }else{
        bool dir = (node->value > pivot->value);

        if(dir){
            curNode->right = node;
            node->parent = curNode;
        }else{
            curNode->left = node;
            node->parent = curNode;
        }

        TreeNode* son = (dir)?pivot->right:pivot->left;
        balanceSingleRotation(pivot, son, dir);

        while(nextNode != ((dir)?pivot->right:pivot->left)){
           nextNode->balance += (node->value > nextNode->value)?1:-1; 
        }
        pivot->balance = 0;
    }
}

void AVLTree::balanceDoubleRotation(TreeNode* pivot, TreeNode* son, TreeNode* grandson, bool dir){
    // rotation on grandson
    if(!dir){
        son->right = grandson;
        grandson->parent = son;

        son->left = grandson->right;
        if(grandson->right)
            grandson->right->parent = son;

        grandson->right = son;
        son->parent = grandson;
    }else{
        son->left = grandson;
        grandson->parent = son;

        son->right = grandson->left;
        if(grandson->left)
            grandson->left->parent = son;

        grandson->left = son;
        son->parent = grandson;

    }

    // rotation on pivot
    if(dir){
        pivot->right = son;
        son->parent = pivot;

        pivot->left = son->right;
        if(son->right)
            son->right->parent = pivot;

        son->right = pivot;
        pivot->parent = son;
    }else{
        pivot->left = son;
        son->parent = pivot;

        pivot->right = son->left;
        if(son->left)
            son->left->parent = pivot;

        son->left = pivot;
        pivot->parent = son;
    }
    // update balances
    pivot->balance = 0; 
}

void AVLTree::balanceSingleRotation(TreeNode* pivot, TreeNode* son, bool dir){
    if(dir){
        pivot->parent->right = son;
        son->parent = pivot->parent;

        pivot->right = son->left;
        if(son->left)
            son->left->parent = pivot;

        son->left = pivot;
        pivot->parent = son;
    }else{
        pivot->parent->left = son;
        son->parent = pivot->parent;

        pivot->left = son->right;
        if(son->right)
            son->right->parent = pivot;

        son->right = pivot;
        pivot->parent = son;
    }
    // update balances
    pivot->balance += (dir)?-1:1;
}

void print(TreeNode* root, int depth = 0){
    print(root->right, depth+1);
    if(root == NULL) return;

    for(int i = 0; i < depth; i++)
        std::cout << "  ";

    std::cout << "->" << root->value << std::endl;
    
    print(root->left, depth+1);
}

MaxHeap::MaxHeap(int capacity){
    heap = new int[capacity];
    this->capacity = capacity;
    size = -1;
}

MaxHeap::~MaxHeap(){
    delete[] heap;
}

void MaxHeap::insert(int value){
    size++;
    if(size >= capacity){
        std::cerr << "size exceeds capacity" << std::endl;
        return;
    }
    heap[size] = value;
    heapify(size);
}

int* MaxHeap::heapsort(){
    int curIdx = size;
    while(curIdx > 0){
        heap[curIdx] ^= heap[0]; 
        heap[0] ^= heap[curIdx]; 
        heap[curIdx] ^= heap[0]; 
        curIdx--;
        int newIdx = 0;
        while(newIdx < curIdx){
            if(heap[curIdx] < heap[2*curIdx + 1] || heap[curIdx] < heap[2*curIdx + 2]){
                if(heap[2*curIdx + 1]  > heap[2*curIdx + 2]){
                    heap[curIdx] ^= heap[2*curIdx + 1]; 
                    heap[2*curIdx + 1] ^= heap[curIdx]; 
                    heap[curIdx] ^= heap[2*curIdx + 1]; 
                    curIdx = 2*curIdx + 1;
                }else{
                    heap[curIdx] ^= heap[2*curIdx + 2]; 
                    heap[2*curIdx + 2] ^= heap[curIdx]; 
                    heap[curIdx] ^= heap[2*curIdx + 2]; 
                    curIdx = 2*curIdx + 2;
                }
            }
        }
    }
    return heap;
}

void MaxHeap::heapify(int from){
    int curIdx = from;
    while(curIdx > 0){
        if(heap[curIdx] > heap[((curIdx - 1)/2)]){
            heap[curIdx] ^= heap[((curIdx - 1)/2)]; 
            heap[((curIdx - 1)/2)] ^= heap[curIdx]; 
            heap[curIdx] ^= heap[((curIdx - 1)/2)]; 
        }
        curIdx = ((curIdx - 1)/2);
    }
}

