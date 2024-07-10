class ListNode:
    def __init__(self, val=0, next=None):
        self.val = val
        self.next = next
    

def rotateRight(head, k: int):        
        print("ptr established")
        length = 0
        ptr1 = head
        ptr2 = head
        jump=head
        while jump!= None:
            length += 1
            jump = jump.next
            
        k = k % length
        position = length - k
        index = 0
        
        if length == 0:
            return head 
        
        while ptr2 != None:
            if index == position:
                break
            else:
                ptr2 = ptr2.next
            index+=1

        
        print("ptr established")
        
        while ptr1.next != None:
            while ptr2 != None:
                temp = ptr2.val
                ptr2.val = ptr1.val
                ptr1.val = temp
                ptr1 = ptr1.next
                ptr2 = ptr2.next
            
            ptr1.next.val, ptr1.val = ptr1.val, ptr1.next.val
            ptr1=ptr1.next
        return head 

myList=ListNode(1)
jump=myList
for i in range(2,6):
    jump.next=ListNode(i)
    jump=jump.next

print(rotateRight(myList,2))



