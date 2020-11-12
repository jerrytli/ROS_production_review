# Recent advancements in ROS compartmentalization, perception, and signaling in plants

Bardo Castro<sup>1</sup>, Matteo Citterico<sup>2</sup>, Sachie Kimura<sup>3</sup>, Danielle M. Stevens<sup>1</sup>, Michael Wraczek<sup>2</sup>, Gitta Coaker<sup>1</sup>.


<sup>1</sup>Department of Plant Pathology, University of California, Davis, USA <br />
<sup>2</sup>Ritsumeikan Global Innovation Research Organization, Ritsumeikan University, Shiga, Japan <br />
<sup>3</sup>Organismal and Evolutionary Biology Research Programme, Viikki Plant Science Centre, Helsinki University, Finland <br />


-----------------------

#### Purpose: The scripts in this repository is for processing and ploting scores of similarity for NADPH oxidase homologs. See the image below or for more details, check out this markdown file [here](process_files.md).


### Part 1: Calculate similarity across full-length NADPH oxidase homologs, the N-terminus, and the C-terminus
Briefly, we wanted to cross compare the similarity of full-length NADPH oxidase homologs to the N- and C-terminus region. This region from NADPH oxidase homologs was extracted via a custom blast search, aligned to their respective sequence comaprsion (full length homolog vs. full length AtRBOHD and parsed N- and C-terminus of homolog vs. N- and C-terminus of AtRBOHD, respectively), and a similarity score was calculated. Unsurprisingly, C-terminus is overall more conserved in similarity to the reference than the N-terminus and full-length protein. 

<img src="https://github.com/DanielleMStevens/ROS_production_review/blob/master/Images/Part1_methods_drawing.png" width="500" height="250">


### Part 2: Calculate similarity across residues of the N- and C-terminus of AtRBOHD
While the above analysis is a simple method to show the C-terminus is more conserved than the rest of the protein, recent work has shown key S/T residues are conserved as part of their role in protein regulation via post-translational modifications. To better assess the similarity within each terminus, a scanning window approach was implimented to assess similarity at a positional basis. Essentially, parsed fragments were locally aligned to the reference and a similarity score was calculated at the reference position.


<img src="https://github.com/DanielleMStevens/ROS_production_review/blob/master/Images/Part2_methods_drawing.png" width="380" height="200">

If you have any questions, concerns, or found a bug, please contact me at dmstev@ucdavis.edu.
