def permute(nums):
        output=[]
        for i in range(len(nums)):
                nums[i]=str(nums[i])
        input=''.join(nums)
        print(input)
        def backTrack(myString,comb):
                
                if(len(comb)==len(nums)):
                        output.append(comb)
                        return
            
                for i in range(len(myString)):
                        backTrack(myString[0:i]+myString[i+1:],comb+myString[i])
                        

            
        backTrack(input,"")
        newOut=[]
        for comb in output:
                temp=[]
                for c in comb:
                        temp.append(int(c))
                newOut.append(temp)
        return newOut


print(permute([1,2,3]))