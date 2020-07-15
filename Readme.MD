# SeqDB Uniprot Local Mirror

Project containing scripts to build and access a local database containing UniprotKB protein info 
via a Python interface using an LMDB or MongoDB backend.

This is used for scripts that need to perform many random queries to fetch protein sequences 
without hitting a web API constantly (slow, probably annoys the server). Essentially this is a 
fancy cache that lets you quickly fetch by multiple types of accessions/IDs.

The LMDB backend is now the default since it is fully embedded and does not require a server.

This now includes a script "seqdb-load" which 