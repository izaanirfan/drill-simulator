def get_annulus_diameter(sections, md):

    for sec in sections:

        if sec.top_md <= md <= sec.end_md:

            # OPEN HOLE
            if sec.type.lower() == "open hole":
                return sec.hole_d

            # CASING
            if sec.type.lower() == "casing":
                return sec.casing_id

            # LINER
            if sec.type.lower() == "liner":

                # below liner top → inside liner
                if md >= sec.top_md:
                    return sec.casing_id

                # above liner → parent casing
                for parent in sections:
                    if parent.type.lower() == "casing" and parent.top_md <= md <= parent.end_md:
                        return parent.casing_id

    return sections[-1].hole_d