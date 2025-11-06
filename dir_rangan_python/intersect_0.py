import torch;
efind = lambda a : torch.where(a)[0] ;

def intersect_0(a_, b_):
    u_a_,index_u_a_from_a_,n_u_a_ = torch.unique(a_.ravel(),return_inverse=True,return_counts=True,dim=0);
    u_b_,index_u_b_from_b_,n_u_b_ = torch.unique(b_.ravel(),return_inverse=True,return_counts=True,dim=0);
    u_a_u_b_ = torch.concatenate([u_a_.ravel(), u_b_.ravel()],0);
    u_u_a_u_b_ , index_u_u_a_u_b_from_u_a_u_b_ , n_u_u_a_u_b_ = torch.unique(u_a_u_b_.ravel(),return_inverse=True,return_counts=True,dim=0);
    index_cap_ = efind(n_u_u_a_u_b_> 1);
    cap_ = u_u_a_u_b_.ravel()[index_cap_];
    n_cap = cap_.numel();
    index_na_from_ncap_ = torch.zeros(n_cap).to(dtype=torch.int32);
    index_nb_from_ncap_ = torch.zeros(n_cap).to(dtype=torch.int32);
    for ncap in range(n_cap):
        cap = cap_[ncap].item();
        na = int(torch.min(efind(a_==cap)).item());
        index_na_from_ncap_[ncap] = na;
        nb = int(torch.min(efind(b_==cap)).item());
        index_nb_from_ncap_[ncap] = nb;
    #end;%for ncap in range(n_cap):
    return(cap_,index_na_from_ncap_,index_nb_from_ncap_);
