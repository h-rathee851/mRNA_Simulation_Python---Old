def __Sim_run(self, float t = 0,np.ndarray pos_jump_cy, np.ndarray pos_w_cy, np.ndarray rib_pos_cy, bool steady_reached = False):
        
        cdef bool enter = False
        cdef bool exit = False
        cdef bool move = False
        
#         cdef int rib_pos_index
        
        
        cdef double tau = Get_tau(self.pos_w_)
        
        t += tau
    
        cdef Py_ssize_t index = Get_mu(self.pos_w_)
        

        
        if steady_reached:
            
            self.__calcDensity(tau)
    
        
        if self.pos_jump_[index] == 0:
            

        
            self.rib_pos_ = np.insert(self.rib_pos_,0,2)

            if len(self.rib_pos_) > 1:
                if self.rib_pos_[1] > self.L+2: 
                    self.pos_jump_[index] = 2
                    self.pos_w_[index] = self.w_[2]
                elif self.rib_pos_[1] <= self.L+2:
                    self.pos_jump_= np.delete(self.pos_jump_, index)
                    self.pos_w_= np.delete(self.pos_w_, index)
            else:
                self.pos_jump_[index] = 2
                self.pos_w_[index] = self.w_[2]
                
            enter = True
        
        elif self.pos_jump_[index] == self.LEN:
        
            self.rib_pos_ = np.delete(self.rib_pos_, -1)
            
            
            if len(self.rib_pos_) >= 1 and self.rib_pos_[-1] == self.LEN-self.L:
                self.pos_jump_[index] = self.rib_pos_[-1]
                self.pos_w_[index] = self.w_[self.rib_pos_[-1]]
            
            else:
                self.pos_jump_= np.delete(self.pos_jump_, index)
                self.pos_w_= np.delete(self.pos_w_, index)
                
            exit = True
        
        else:
        
            rib_pos_index = np.where(self.rib_pos_ == self.pos_jump_[index])
            rib_pos_index = rib_pos_index[0][0]
            self.rib_pos_[rib_pos_index] += int(1)
            
            
            if len(self.rib_pos_) > rib_pos_index+1:
                if self.rib_pos_[rib_pos_index+1] - self.rib_pos_[rib_pos_index] > self.L:
                    self.pos_jump_[index] = self.rib_pos_[rib_pos_index]
                    self.pos_w_[index] = self.w_[self.rib_pos_[rib_pos_index]]
                else:
                    self.pos_jump_= np.delete(self.pos_jump_, index)
                    self.pos_w_= np.delete(self.pos_w_, index)
            else:
                
                self.pos_jump_[index] = self.rib_pos_[rib_pos_index]
                self.pos_w_[index] = self.w_[self.rib_pos_[rib_pos_index]]
                
            if rib_pos_index != 0 and self.rib_pos_[rib_pos_index] - self.rib_pos_[rib_pos_index-1] == self.L+1 :
                
                self.pos_jump_ = np.insert(self.pos_jump_,index,self.rib_pos_[rib_pos_index-1]) 
                self.pos_w_ = np.insert(self.pos_w_,index,self.w_[self.rib_pos_[rib_pos_index-1]])
            
                
            elif rib_pos_index == 0 and self.rib_pos_[rib_pos_index] == self.L+2:
                self.pos_jump_ = np.insert(self.pos_jump_,0,0) 
                self.pos_w_ = np.insert(self.pos_w_,0,self.w_[0])
        
            move = True
    
 
        return t, enter, exit, move, pos_jump_cy, pos_w_cy, rib_pos_cy 
