function json_str = get_or_pitch_roll(or_list_filepath)

if(exist(or_list_filepath, 'file') ~= 2)
    error('404: file not found')
end

orl_data = ORReader(or_list_filepath);
ors = orl_data.ORs;

sched_start = time([2021, 01, 25, 00, 00, 00]);
sched_stop = sched_start + 86400*7;
delta_t = 3280;
sched_t = sched_start:delta_t:sched_stop;

eph = EphemerisEngine('target', 'Sun');
eph.time= sched_t;
sun_pos=eph.Position;
characteristics = axaf_characteristics;

pitch_deg = cell(1, numel(ors));
roll_deg = pitch_deg;
dwell_t = NaN(1, numel(ors));
or_simpos = NaN(1, numel(ors));
ccd_cnt = NaN(1, numel(ors));
opt_chips = NaN(1, numel(ors));

for idx = 1:numel(ors)
    %Find a target attitude for this OR.
    dwell_t(idx) = ors(idx).duration(1);
    ra = ors(idx).target(1)*pi/180;
    dec = ors(idx).target(2).*pi/180;
    offY = ors(idx).targetOffset(1)*pi/180;
    offZ = ors(idx).targetOffset(2)*pi/180;
    si_field = strrep(ors(idx).SI,'-','_');
    M = characteristics.Science.ACA.SI2ACA.(si_field);
    or_simpos(idx) = characteristics.Science.(si_field).simpos;
    ccd_cnt(idx) = ors(idx).ReqChip;
    opt_chips(idx) = ors(idx).OptChip;
    
    %First find a nominal roll attitude for this target
    [att, nominal_roll] = aca_matrix(ra, dec, offY, offZ, M, 0, sun_pos); %#ok<ASGLU>
%     % simplified calculation not used here
%     radr = [repmat([ra, dec], numel(nominal_roll), 1), nominal_roll'];
%     att = NaN(3,3,numel(nominal_roll));
%     for roll_idx = 1:numel(nominal_roll)
%         att(:, :, roll_idx) = radr2att(radr(roll_idx, :));
%     end
    quats = att2quat(att);
    
    [pitch, roll] = pitchroll(q_rotate(quats, unit_vec(sun_pos)));
    
    pitch_deg{idx} = pitch .* 180/pi;
    roll_deg{idx} = roll .* 180/pi;
    non_nominal_roll = any(abs(roll_deg{idx}) > 1e-10);
    if(non_nominal_roll)
        error('calculated roll was farther off-nominal than expected.')
    end
end

to_py.pitch_deg = pitch_deg;
to_py.sched_t = cellstr(sched_t);
to_py.dwell_t = dwell_t;
to_py.or_simpos = or_simpos;
to_py.ccd_cnt = ccd_cnt;
to_py.opt_chips = opt_chips;
json_str = jsonencode(to_py);
end
