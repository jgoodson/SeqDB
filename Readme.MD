Project containing scripts to build and access a local database containing UniprotKB protein info via a Python interface using a MongoDB backend.

This is used for scripts that need to perform many random queries to fetch protein sequences without hitting a web API constantly (slow, probably annoys the server). Essentially this is a fancy cache that lets you quickly fetch by multiple types of accessions/IDs.

This version also includdes a SeqVecDB module containing code for a parallel MongoDB containing vector representations of protein sequences embedding using one or more method.

Currently, this includes code to download and utilize a pre-trained model from https://github.com/stamakro/GCN-for-Structure-and-Function