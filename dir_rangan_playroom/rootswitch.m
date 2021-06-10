function string_out_ = rootswitch(string_0in_,string_root_,string_key_);
string_out_ = sprintf('/%s/%s',string_root_,string_0in_(strfind(string_0in_,string_key_):end));
