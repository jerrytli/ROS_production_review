# Recent advancements in ROS production, regulation, and signaling in plants

Bardo Castro<sup>1</sup>, Matteo Citterico<sup>2</sup>, Sachie Kimura<sup>3</sup>, Danielle M. Stevens<sup>1</sup>, Michael Wraczek<sup>2</sup>, Gitta Coaker<sup>1</sup>.


<sup>1</sup>Department of Plant Pathology, University of California, Davis, USA <br />
<sup>2</sup>Ritsumeikan Global Innovation Research Organization, Ritsumeikan University, Shiga, Japan <br />
<sup>3</sup>Organismal and Evolutionary Biology Research Programme, Viikki Plant Science Centre, Helsinki University, Finland <br />


-----------------------

Purpose: The script in this repository is for processing and ploting scores of similarity for NADPH oxidase homologs.


### Part 1: Calculate global full-length NADPH oxidase homologs 
Briefly, to cross compare the similarity of full-length NADPH oxidase homologs to the more conserved C-terminus region. This region from NADPH oxidase homologs, aligned to their respective sequence comaprsion (full length homolog vs. full length AtRBOHD and parse C-terminus of homolog vs. C-terminus of AtRBOHD). Unsurprisingly, this region is overall more convered in similarity to the reference than the full-length protein. See the image below or for more details, check out this markdown file [here]().

### Part 2: Calculate similarity across residues of C-terminus AtRBOHD
While the above analysis is a simple method to show the C-terminus is more conserved than the rest of the protein, recent work has shown key S/T residues are conserved as part of their role in post-translational control. To better assess the similarity within this terminus, a scanning window approach was implimented with a 20 amino acid window to assess similarity at a positional basis. See the image below or for more details, check out this markdown file [here]().
