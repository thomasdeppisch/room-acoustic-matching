---
layout: default
---

<script src="https://cdn.rawgit.com/download/polymer-cdn/1.5.0/lib/webcomponentsjs/webcomponents-lite.min.js"></script>
<link href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet" integrity="sha384-wvfXpqpZZVQGK6TAh5PVlGOfQNHSoD2xbE+QkPxCAFlNEevoEH3Sl0sibVcOQVnN" crossorigin="anonymous">
<link rel="stylesheet" href="thirdparty/trackswitch-js/trackswitch.min.css" />
<script src="https://code.jquery.com/jquery-3.2.1.min.js" integrity="sha256-hwg4gsxgFZhOsEEamdOYGBf13FyQuiTwlAQgxVSNgt4=" crossorigin="anonymous"></script>
<script src="thirdparty/trackswitch-js/trackswitch.min.js"></script>
<script type="text/javascript">
    var settings = {
        onlyradiosolo: true,
        repeat: true,
    };

    jQuery(document).ready(function() {
        jQuery(".player").trackSwitch(settings); 
    });
</script>

# Identification and Matching of Room Acoustics With Moving Head-Worn Microphone Arrays
This page presents binaural audio examples for the manuscript
> T. Deppisch, S. Amengual Garí, P. Calamia, and J. Ahrens, "Identification and Matching of Room Acoustics With Moving Head-Worn Microphone Arrays", 2026.

## Abstract
Head-worn devices such as smartglasses and headsets are the predominant form factor for augmented reality and telepresence applications, where real-world environments are augmented with virtual sound sources. For these sources to appear perceptually convincing, the acoustics of their virtual environment must closely match the room acoustics of the physical space. Estimating room impulse responses (RIRs) in this setting is challenging because practical scenarios require small arrays that continuously move with the user's head. This study presents a method for the blind identification of RIRs from speech signals captured with a moving head-worn microphone array and the subsequent rendering of virtual sound sources based on perceptually relevant acoustic parameter estimates. A motion-aware signal model that estimates spatial RIRs as sound field coefficients and incorporates position tracking data is compared against an omnidirectional model and a baseline. Numerical results show that the motion-aware model provides the most accurate acoustic parameter estimates when used in an informed setting with the true reference signal. In the blind setting, however, its advantage largely diminishes. In a listening experiment, renderings based on the omnidirectional model are rated as most similar to the reference condition and are most often associated with the correct room significantly above chance. The findings highlight the practical relevance of the proposed framework, with the omnidirectional model offering robust and perceptually convincing performance, while the motion-aware model remains promising for parameter estimation.

## Examples
Below are examples from both parts of the listening experiment described in the paper.

### MUSHRA
{% assign example_ids = "1,2,3,4,5,6,7,8" | split: "," %}
{% for example_id in example_ids %}
### Example {{ forloop.index }}

<div class="player">
    <ts-track title="Reference">
        <ts-source src="audio_examples/1_MUSHRA/{{ example_id }}_ref.wav"></ts-source>
    </ts-track>
    <ts-track title="Hidden Reference">
        <ts-source src="audio_examples/1_MUSHRA/{{ example_id }}_hiddenRef.wav"></ts-source>
    </ts-track>
    <ts-track title="Anchor">
        <ts-source src="audio_examples/1_MUSHRA/{{ example_id }}_anchor.wav"></ts-source>
    </ts-track>
    <ts-track title="Estimated Full Model">
        <ts-source src="audio_examples/1_MUSHRA/{{ example_id }}_estFull.wav"></ts-source>
    </ts-track>
    <ts-track title="Estimated Non-Blind Model">
        <ts-source src="audio_examples/1_MUSHRA/{{ example_id }}_estNonBlind.wav"></ts-source>
    </ts-track>
    <ts-track title="Estimated Simple Model">
        <ts-source src="audio_examples/1_MUSHRA/{{ example_id }}_estSimple.wav"></ts-source>
    </ts-track>
</div>

{% endfor %}

### 2AFC
{% assign afc_examples = "1|office|hiddenRef|Hidden Reference|storage,2|office|hiddenRef|Hidden Reference|kitchen,6|storage|hiddenRef|Hidden Reference|hallway,7|kitchen|hiddenRef|Hidden Reference|office,12|hallway|hiddenRef|Hidden Reference|kitchen,25|office|est|Estimated|storage,30|storage|est|Estimated|hallway,31|kitchen|est|Estimated|office,35|hallway|est|Estimated|storage,46|hallway|est|Estimated|office" | split: "," %}
{% for example in afc_examples %}
{% assign fields = example | split: "|" %}
{% assign trial_id = fields[0] %}
{% assign reference_room = fields[1] %}
{% assign test_condition = fields[2] %}
{% assign test_title = fields[3] %}
{% assign other_room = fields[4] %}
### Example {{ forloop.index }}: {{ reference_room | capitalize }} vs. {{ other_room | capitalize }}

<div class="player">
    <ts-track title="Reference ({{ reference_room | capitalize }})">
        <ts-source src="audio_examples/2_2AFC/{{ trial_id }}a_ref_{{ reference_room }}.wav"></ts-source>
    </ts-track>
    <ts-track title="{{ test_title }} ({{ reference_room | capitalize }})">
        <ts-source src="audio_examples/2_2AFC/{{ trial_id }}b_{{ test_condition }}_{{ reference_room }}.wav"></ts-source>
    </ts-track>
    <ts-track title="Other Room ({{ other_room | capitalize }})">
        <ts-source src="audio_examples/2_2AFC/{{ trial_id }}c_other_{{ other_room }}.wav"></ts-source>
    </ts-track>
</div>

{% endfor %}
